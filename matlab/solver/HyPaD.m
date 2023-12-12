function [L,U,N,ids,it,exitflag] = HyPaD(problem,param,HYPAD_PARAM,L,U,EPSILON,OFFSET)

%% Initialization phase
% General initialization
[n,m,p,q,f,g,Df,Dg,Aineq,bineq,Aeq,beq,lb,ub,x0,~,~] = problem(param);
N = [];
E = [];
ids = cell(0,5);
it = 0;
exitflag = -3;
TIME_LIMIT = 3600;

% SNIA Initialization
ALL_INTEGERS_VISITED = false;
SNIA_MODE = HYPAD_PARAM(1);
SNIA_OPTIONS = HYPAD_PARAM(2);
SNIA_GUESS = HYPAD_PARAM(3);
if SNIA_MODE == 1
    [boxesLB,boxesUB,int_num,current_int]=SNIABoxesInit(n,m,lb,ub,SNIA_OPTIONS);
elseif SNIA_MODE == 2
    [boxesLB,boxesUB,int_num,current_int]=SNIABoxesInit(n,m,lb,ub,0);
elseif SNIA_MODE == 3
    [int_lb,width,int_num,current_int]=SNIAListInit(n,m,lb,ub);
end

%% Check the initial solution
sol_x_integer = round(x0(n+1:n+m));
[~,err] = FeasibilityCheckPatch(n,q,g,Aineq,bineq,Aeq,beq,lb,ub,x0,sol_x_integer);
if err<1e-6
    [eff_candidates,nd_candidates] = WSNLP(n,m,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,sol_x_integer,x0,eye(p));
    z = min(nd_candidates,[],2);
    N = [N nd_candidates];
    E = [E eff_candidates];
    for k=1:size(nd_candidates,2)
        U = updateLUB3(U,nd_candidates(:,k));
    end
    x0 = eff_candidates(:,ceil(p/2));
    ids(end+1,:) = {sol_x_integer,z-OFFSET,true,nd_candidates,x0};
    X = [];
else
    x0 = [x0(1:n);sol_x_integer];
    X = [x0];
end

%% Main Loop
tic;
while (toc<TIME_LIMIT)
    it = it+1;
    size_L = size(L,2);
    finished_bounds = 0;
    rsup_error = true;
    for l=L
        if (toc > TIME_LIMIT), break; end %recognize time limit early
        U_temp = all(l<(U-EPSILON+OFFSET));
        if any(U_temp)
            directions = U(:,U_temp)-l;
%             [~,u_index] = max(prod(directions)); %Hypervolume related
            [~,u_index] = max(min(directions)); %Width related
            r = directions(:,u_index);
            u = l+r;
            [solutions,flag_rsup]=RSUP(n,m,p,q,f,g,Df,Dg,Aineq,bineq,Aeq,beq,lb,ub,x0,[X E],l,r);
            if flag_rsup ~= 1
                warning('HyPaD:RSUP_flag','RSUP terminated with exitflag ~= 1');
                if flag_rsup == -3
                    exitflag = -1;
                    break;
                end
            else
                sol_x_integer = round(solutions(n+1:n+m));
                sol_eta = solutions(n+m+2:n+m+1+p);
                rsup_error = false;
                if sol_eta > u+OFFSET
                    warning('HyPaD:RSUP_eta','RSUP returned eta > u');
                    [sol_relax,~,flag_rsup,~] = PS(f,g,Aineq,bineq,Aeq,beq,lb,ub,x0,l,r);
                    if flag_rsup < 0
                        warning('HyPaD:RSUPIR_flag','RSUPIR returned infeasible solution');
                    else
                        sol_relax_eta = f(sol_relax);
                        L = updateLLB3(L,sol_relax_eta);
                    end
                else
                    L = updateLLB3(L,sol_eta);
                end
                if ~ALL_INTEGERS_VISITED
                    [x,err] = FeasibilityCheckPatch(n,q,g,Aineq,bineq,Aeq,beq,lb,ub,x0,sol_x_integer);
                    if err < OFFSET
                        if ~isempty(ids)
                            if m > 1
                                ids_index = find(all([ids{:,1}]==sol_x_integer));
                            else
                                ids_index = find([ids{:,1}]==sol_x_integer);
                            end
                        else
                            ids_index = [];
                        end
                        if isempty(ids_index)
                            %% InitIDS
                            [eff_candidates,nd_candidates] = WSNLP(n,m,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,sol_x_integer,x0,eye(p));
                            z = min(nd_candidates,[],2);
                            N = [N nd_candidates];
                            E = [E eff_candidates];
                            for k=1:size(nd_candidates,2)
                                U = updateLUB3(U,nd_candidates(:,k));
                            end
                            x0 = eff_candidates(:,ceil(p/2));
                            ids(end+1,:) = {sol_x_integer,z-OFFSET,true,nd_candidates,x0};
                        else
                            if ~ids{ids_index,3} %selected patch inactive
                                ids_index = find([ids{:,3}],1);
                            end
                            if ~isempty(ids_index)
                                %% UpdateIDS
                                patch_L = ids{ids_index,2};
                                patch_done = true;
                                x0 = ids{ids_index,5};
                                for patch_l=patch_L
                                    if (toc > TIME_LIMIT), break; end %recognize time limit early
                                    patch_U_temp = all(patch_l<(U-EPSILON+OFFSET));
                                    if any(patch_U_temp)
                                        patch_done = false;
                                        patch_directions = U(:,patch_U_temp)-patch_l;
                                        [~,patch_u_index] = max(min(patch_directions));
                                        d = patch_directions(:,patch_u_index);
                                        [x_sup,t_sup,flag_sup,~] = SUP(n,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,ids{ids_index,1},x0,patch_l,d);
                                        if flag_sup < 0
                                            warning('HyPaD:SUP_flag','Infeasible solution obtained for SUP');
                                        else
                                            y_lub = f(x_sup);
                                            y_llb = patch_l+t_sup.*d;
                                            ids{ids_index,2} = updateLLB3(ids{ids_index,2},y_llb);
                                            U = updateLUB3(U,y_lub);
                                            ids{ids_index,4} = [ids{ids_index,4}, y_lub];
                                            N = [N y_lub];
                                            E = [E x_sup];
                                            x0 = x_sup;
                                        end 
                                    end
                                end
                                ids{ids_index,5} = x0;
                                if patch_done
                                    ids{ids_index,3} = false;
                                end
                            else
                                %% SNIA
                                if SNIA_MODE == 1
                                    [x,SNIA_flag,ALL_INTEGERS_VISITED,current_int] = SNIABoxesFixed(n,m,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,ids,x0,X,boxesLB,boxesUB,int_num,current_int,l,r,OFFSET,SNIA_GUESS);
                                elseif SNIA_MODE == 2
                                    [x,SNIA_flag,ALL_INTEGERS_VISITED,current_int] = SNIABoxesDynamic(n,m,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,ids,x0,X,boxesLB,boxesUB,int_num,l,r,OFFSET,SNIA_OPTIONS,SNIA_GUESS);
                                elseif SNIA_MODE == 3
                                    [x,SNIA_flag,ALL_INTEGERS_VISITED,current_int] = SNIAList(n,m,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,ids,x0,X,int_lb,width,int_num,current_int,l,r,OFFSET,SNIA_GUESS);
                                end
                                if SNIA_flag > 0
                                    %% InitIDS
                                    sol_x_integer = x(n+1:n+m);
                                    [eff_candidates,nd_candidates] = WSNLP(n,m,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,sol_x_integer,x0,eye(p));
                                    z = min(nd_candidates,[],2);
                                    N = [N nd_candidates];
                                    E = [E eff_candidates];
                                    for k=1:size(nd_candidates,2)
                                        U = updateLUB3(U,nd_candidates(:,k));
                                    end
                                    x0 = eff_candidates(:,ceil(p/2));
                                    ids(end+1,:) = {sol_x_integer,z-OFFSET,true,nd_candidates,x0};
%                                     disp('INFO: SNIA computed a new feasible integer assignment.');
                                else
                                    if ALL_INTEGERS_VISITED
                                        %% Terminate HyPaD
                                        disp('END:  All integer assigments visited! Switching solver mode.');
                                        L = [ids{:,2}];
                                        L_indexlist = computeNDS(L);
                                        L = L(:,L_indexlist);
                                        exitflag = 0;
                                        return;
                                    else
                                        if SNIA_flag < 0
%                                             disp('INFO: SNIA computed an invalid integer assignment.');
                                        else
                                            x0 = x;
                                            X = [X x];
%                                             disp('INFO: SNIA computed a new infeasible integer assignment.');
                                        end
                                    end
                                end
                            end
                        end
                    else
                        x0 = x;
                        X = [X x];
                    end
                end
            end
            indexlist = [ids{:,3}];
        else
            finished_bounds = finished_bounds+1;
        end
    end
    if finished_bounds > size_L-0.5
        disp("END:  Finished! Width is <= EPSILON.");
        exitflag = 1;
        return;
    end
    if rsup_error
        disp("ERR:  Algorithm terminated due to RSUP issues.");
        exitflag = -2;
        return;
    end
    if ~any(indexlist)
%         disp('INFO: Iteration ended with all patches inactive');
    end
end
if toc > TIME_LIMIT
    exitflag = -1;
end
end