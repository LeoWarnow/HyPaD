function [L,U,N,ids,it,exitflag,time] = callSolver(problem_name,problem_param,L,U,HYPAD_PARAM,EPSILON,OFFSET,plot_result)

%% Load problem
problem = str2func(problem_name);

%% Handling of incorrect or missing input
[~,m,p,~,f,~,~,~,~,~,~,~,lb,ub,~,is_convex,is_quadratic] = problem(problem_param);

% Applying standard tolerances
if isempty(EPSILON)
    EPSILON = 0.1;
end
if isempty(OFFSET)
    OFFSET = EPSILON*1e-3;
end

% Checking HyPaD parameters
if isempty(HYPAD_PARAM)
    HYPAD_PARAM = [1,4,0];
    disp('Selected SNIA with fixed boxes, 4 splits and no guess.');
end

% Initialization of L,U
if isempty(L) || isempty(U)
    startbox = infsup(lb, ub);
    B = intval;
    for i=1:p
        B(i) = (1:p==i)*f(startbox);
    end
    if isempty(L)
        L = B.inf';
    end
    if isempty(U)
        U = B.sup';
    end
end
L = L-OFFSET;
U = U+OFFSET;

%% Call HyPaD
tic;
[L,U,N,ids,it,exitflag] = HyPaD(problem,problem_param,HYPAD_PARAM,L,U,EPSILON,OFFSET);
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