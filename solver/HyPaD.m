function [L,U,N,ids,it,hypad_flag] = HyPaD(problem,param,x_init,z,Z,EPSILON,OFFSET)

%% Initialization phase
% General initialization
[n,m,p,q,f,g,Df,Dg,Aineq,bineq,Aeq,beq,lb,ub] = problem(param);
N = [];
E = [];
B = [];
x0 = x_init;
L = z;
U = Z;
ids = cell(0,5);
it = 0;
hypad_flag = -3;
TIME_LIMIT = 3600;

% SNIA Initialization
ALL_INTEGERS_VISITED = false;
[boxesLB,boxesUB,int_num,current_int]=initSNIA(n,m,lb,ub,4);

%% Check the initial solution
sol_x_integer = round(x_init(n+1:n+m));
[~,err] = FeasibilityCheck(n,q,g,Aineq,bineq,Aeq,beq,lb,ub,sol_x_integer,x0);
if err<1e-6
    [eff_candidates,nd_candidates] = WSNLP(n,m,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,sol_x_integer,x0,eye(p));
    [z,~]=initBoundsNLP(n,m,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,sol_x_integer,x0);
    N = [N nd_candidates];
    E = [E eff_candidates];
    for k=1:size(nd_candidates,2)
        U = updateLUB3(U,nd_candidates(:,k));
    end
    x0 = eff_candidates(:,ceil(p/2));
    ids(end+1,:) = {sol_x_integer,z-OFFSET,true,nd_candidates,x0};
    X = [];
else
    x0 = [x_init(1:n);sol_x_integer];
    X = [x0];
end

%% Main Loop
tic;
while (toc<TIME_LIMIT)
    it = it+1;
    size_L = size(L,2);
    size_U = size(U,2);
    L_temp = L;
    finished_bounds = 0;
    error_counter = 0;
    for l=L_temp
        U_temp = all((repmat(l,1,size_U)<U-EPSILON));
        if any(U_temp)
            directions = U(:,U_temp)-l;
                [~,u_index] = max(min(directions));
                r = directions(:,u_index);
                u = l+r;
                [solutions,flag]=RSUP(n,m,p,q,f,g,Df,Dg,Aineq,bineq,Aeq,beq,lb,ub,[X E],l,r);
                if flag ~= 1
                    warning('HyPaD:RSUP_flag','RSUP terminated with exitflag ~= 1!');
                    error_counter = error_counter+1;
                else
                    sol_x_integer = round(solutions(n+1:n+m));
                    sol_eta = solutions(n+m+2:n+m+1+p);
                    if sol_eta >= u
                        warning('HyPaD:RSUP_eta','RSUP returned eta > u!');
                        sol_relax = RSUPIR(n,m,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,x0,l,r);
                        sol_relax_eta = f(sol_relax);
                        L = updateLLB3(L,sol_relax_eta);
                    else
                        L = updateLLB3(L,sol_eta);
                    end
                    if ~ALL_INTEGERS_VISITED
                        [x,err] = FeasibilityCheck(n,q,g,Aineq,bineq,Aeq,beq,lb,ub,sol_x_integer,x0);
                        if err < 1e-6
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
                                [z,~]=initBoundsNLP(n,m,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,sol_x_integer,x0);
                                N = [N nd_candidates];
                                E = [E eff_candidates];
                                for k=1:size(nd_candidates,2)
                                    U = updateLUB3(U,nd_candidates(:,k));
                                end
                                size_U = size(U,2);
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
                                        patch_U_temp = all((repmat(patch_l,1,size_U)<U-EPSILON));
                                        if any(patch_U_temp)
                                            patch_done = false;
                                            patch_directions = U(:,patch_U_temp)-patch_l;
                                            [~,patch_u_index] = max(min(patch_directions));
                                            d = patch_directions(:,patch_u_index);
                                            [t,x_sup] = SUP(n,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,ids{ids_index,1},x0,patch_l,d); %parfeval?
                                            y_lub = f(x_sup);
                                            y_llb = patch_l+t.*d;
                                            if ~(y_llb < Z)
                                                t = min((Z-patch_l)./d);
                                                y_llb = patch_l+t.*d-OFFSET;
                                                warning('HyPaD:Patch','Update point outside of box [z,Z]!');
                                            end
                                            ids{ids_index,2} = updateLLB3(ids{ids_index,2},y_llb);
                                            U = updateLUB3(U,y_lub);
                                            size_U = size(U,2);
                                            ids{ids_index,4} = [ids{ids_index,4}, y_lub];
                                            N = [N y_lub];
                                            E = [E x_sup];
                                            x0 = x_sup;
                                        end
                                    end
                                    ids{ids_index,5} = x0;
                                    if patch_done
                                        ids{ids_index,3} = false;
                                    end
                                else
                                    %% SNIA
                                    [x,B_new,SNIA_flag,ALL_INTEGERS_VISITED,current_int]=SNIA(n,m,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,ids,x0,X,B,boxesLB,boxesUB,int_num,current_int,l,r); %Boxes
                                    if SNIA_flag > 0
                                        %% InitIDS
                                        sol_x_integer = x(n+1:n+m);
                                        [eff_candidates,nd_candidates] = WSNLP(n,m,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,sol_x_integer,x0,eye(p));
                                        [z,~]=initBoundsNLP(n,m,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,sol_x_integer,x0);
                                        N = [N nd_candidates];
                                        E = [E eff_candidates];
                                        for k=1:size(nd_candidates,2)
                                            U = updateLUB3(U,nd_candidates(:,k));
                                        end
                                        size_U = size(U,2);
                                        x0 = eff_candidates(:,ceil(p/2));
                                        ids(end+1,:) = {sol_x_integer,z-OFFSET,true,nd_candidates,x0};
                                        disp('INFO: SNIA computed a new feasible integer assignment.');
                                    else
                                        if ALL_INTEGERS_VISITED
                                            %% Terminate HyPaD
                                            disp('END:  All integer assigments visited! Switching solver mode.');
                                            TempListLLB = [ids{:,2}];
                                            sizeTempLLB = size(TempListLLB,2);
                                            LLB_indexlist = true(1,sizeTempLLB);
                                            for k=1:sizeTempLLB
                                                LLB_indexlist(k) = ~any(all(bsxfun(@ge,TempListLLB(:,k),TempListLLB(:,[1:k-1,k+1:end]))));
                                            end
                                            L = TempListLLB(:,LLB_indexlist);
                                            hypad_flag = 0;
                                            return;
                                        else
                                            if SNIA_flag < 0
                                                disp('INFO: SNIA computed an invalid integer assignment.');
                                            else
                                                x0 = x;
                                                X = [X x];
                                                disp('INFO: SNIA computed a new infeasible integer assignment.');
                                            end
                                        end
                                    end
                                    if ~isempty(B_new)
                                        B = [B, B_new];
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
        hypad_flag = 1;
        return;
    end
    if error_counter > size_L-0.5
        disp("ERR:  Algorithm terminated due to RSUP issues.");
        hypad_flag = -2;
        return;
    end
    if ~any(indexlist)
        disp('INFO: Iteration ended with all patches inactive!');
    end
end
if toc > TIME_LIMIT
    hypad_flag = -1;
end
end