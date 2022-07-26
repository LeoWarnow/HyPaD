function [x,err]=FeasibilityCheck(n,q,g,Aineq,bineq,Aeq,beq,lb,ub,x_int,x_start)
%FeasibilityCheck Checkes if there exists x such that (x,x_int) is feasible

% Initialization
x0 = [x_start(1:n);0];

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
% nonlcon_rhs = zeros(q,1);

options = optimoptions('fmincon','Display','none');
% options = optiset('solver','ipopt','display','off');

solution = fmincon(@(x) x(n+1),x0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
% opti_object = opti('fun',@(x) x(n+1),'ineq',Aineq,bineq,'eq',Aeq,beq,'nlmix',nonlcon,nonlcon_rhs,-ones(q,1),'bounds',lb,ub,'options',options);
% [solution,~,~,~] = solve(opti_object,x0);
x = [solution(1:n);x_int];
err = solution(n+1);
end