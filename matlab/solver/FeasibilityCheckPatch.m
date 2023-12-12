function [x,err]=FeasibilityCheckPatch(n,q,g,Aineq,bineq,Aeq,beq,lb,ub,x_start,x_int)
%FeasibilityCheckPatch Checkes if there exists x such that (x,x_int) is feasible

% If there are any constraints for the original problem then the continuous
% part has to be copied for the subproblem
if ~isempty(Aineq)
    size_Aineq = size(Aineq,1);
    bineq = bineq-Aineq(:,n+1:end)*x_int;
    Aineq = [Aineq(:,1:n),-ones(size_Aineq,1)];
end
if ~isempty(Aeq)
    size_Aeq = size(Aeq,1);
    beq = beq-Aeq(:,n+1:end)*x_int;
    Aeq = [Aeq(:,1:n),zeros(size_Aeq,1)];
end

if ~isempty(lb)
    lb = [lb(1:n);-Inf];
else
    lb = [-Inf(n,1);-1];
end
if ~isempty(ub)
    ub = [ub(1:n);Inf];
end

function [c,ceq] = nonlcon_fun(x)
    c = g([x(1:n);x_int]) - x(n+1).*ones(q,1);
    ceq = [];
end
nonlcon = @nonlcon_fun;
options = optimoptions('fmincon','Display','none');

x0 = [x_start(1:n);0];

solution = fmincon(@(x) x(n+1),x0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
x = [solution(1:n);x_int];
err = solution(n+1);
end