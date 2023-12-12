function [solution_t,solution_x] = SUP(n,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,x_int,x_start,l,d)
%SUP Search Update Point starting from l along a specified direction d
%which is typically (u-l)

% Initialization of parameters for the optimization problem / solver
x0 = [x_start(1:n);0];

% If there are any constraints for the original problem then the continuous
% part has to be copied for the subproblem
if ~isempty(Aineq)
    size_Aineq = size(Aineq,1);
    bineq = bineq-Aineq(:,n+1:end)*x_int;
    Aineq = [Aineq(:,1:n),zeros(size_Aineq,1)];
end
if ~isempty(Aeq)
    size_Aeq = size(Aeq,1);
    beq = beq-Aeq(:,n+1:end)*x_int;
    Aeq = [Aeq(:,1:n),zeros(size_Aeq,1)];
end
if ~isempty(lb)
    lb = [lb(1:n);-Inf];
end
if ~isempty(ub)
    ub = [ub(1:n);Inf];
end
function [c,ceq] = nonlcon_fun(x)
    x_part = [x(1:n);x_int];
    c_new = -l-x(n+1).*d+f(x_part);
    c = [g(x_part);c_new];
    ceq = [];
end
nonlcon_sup = @nonlcon_fun;
options = optimoptions('fmincon','Display','none');

solution = fmincon(@(x) x(n+1),x0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon_sup,options);
solution_t = solution(n+1);
solution_x = [solution(1:n);x_int];
end