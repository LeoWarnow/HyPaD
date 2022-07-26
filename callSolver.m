function [L,U,N,ids,it,flag,time,z,Z] = callSolver(problem_name,param,z,Z,EPSILON,OFFSET,plot_result)

%% Load additional paths and problem
addpath(genpath('problems'));
addpath(genpath('solver'));
problem = str2func(problem_name);
[n,m,p,q,f,g,~,~,Aineq,bineq,Aeq,beq,lb,ub,x0] = problem(param);

%% Initialization of L,U
if isempty(z) || isempty(Z)
    startbox = infsup(lb, ub);
    B = intval;
    for i=1:p
        B(i) = (1:p==i)*f(startbox);
    end
    if isempty(z)
        z = B.inf';
    end
    if isempty(Z)
       Z = B.sup';
    end
end
z = z-OFFSET;
Z = Z+OFFSET;
x_init = RSUPIR(n,m,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,x0,z,Z-z);

%% Call HyPaD
tic;
[L,U,N,ids,it,flag] = HyPaD(problem,param,x_init,z,Z,EPSILON,OFFSET);
time = toc;

%% Plot if wanted
if plot_result > 0
    plotBoxes(L,U,p);
    if p < 3
        figure;
        hold on;
        plot(L(1,:),L(2,:),'LineStyle','none','Marker','.','Color','blue');
        plot(U(1,:),U(2,:),'LineStyle','none','Marker','.','Color',[1, 0.4745, 0]);
        grid on;
        xlabel('f_1');
        ylabel('f_2');
    elseif p < 4
        figure;
        hold on;
        plot3(L(1,:),L(2,:),L(3,:),'LineStyle','none','Marker','.','Color','blue');
        plot3(U(1,:),U(2,:),U(3,:),'LineStyle','none','Marker','.','Color',[1, 0.4745, 0]);
        grid on;
        xlabel('f_1');
        ylabel('f_2');
        zlabel('f_3');
    end
end
end