function [sol_x]=RSUPIR(n,m,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,x_start,a,r)
%RSUPIR Computes a solution of the integer relaxed version of RSUP

% Initialization
x0 = [x_start;0];

% Adjustment of existing constraints
if ~isempty(Aineq)
    sizeAineq = size(Aineq,1);
    Aineq = [Aineq,zeros(sizeAineq,1)];
end
if ~isempty(Aeq)
    size_Aeq = size(Aeq,1);
    Aeq = [Aeq,zeros(size_Aeq,1)];
end
if ~isempty(lb)
    lb = [lb;-Inf];
end
if ~isempty(ub)
    ub = [ub;Inf];
end

function [c,ceq] = nonlcon_fun(x)
    x_part = x(1:n+m);
    c_new = -a-x(n+m+1).*r+f(x_part);
    c = [g(x_part);c_new];
    ceq = [];
end
nonlcon = @nonlcon_fun;
options = optimoptions('fmincon','Display','none');

sol = fmincon(@(x) x(n+m+1),x0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
sol_x = sol(1:n+m);
end