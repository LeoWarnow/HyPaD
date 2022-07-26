function [sol_x,sol_f]=WSNLP(n,m,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,x_int,x_start,weights)
%WSNLP Solves weighted sum scalarizations of NLP

% Initialization of parameters for the optimization problem / solver
x0 = x_start(1:n);

% If there are any constraints for the original problem then the continuous
% part has to be copied for the subproblem
if ~isempty(Aineq)
    bineq = bineq-Aineq(:,n+1:end)*x_int;
    Aineq = Aineq(:,1:n);
end
if ~isempty(Aeq)
    beq = beq-Aeq(:,n+1:end)*x_int;
    Aeq = Aeq(:,1:n);
end
if ~isempty(lb)
    lb = lb(1:n);
end
if ~isempty(ub)
    ub = ub(1:n);
end
function [c,ceq] = nonlcon_fun(x)
    c = g([x;x_int]);
    ceq = [];
end
nonlcon = @nonlcon_fun;
% nonlcon_rhs = zeros(q,1);

% Solver options for OPTI toolbox
options = optimoptions('fmincon','Display','none');
% options = optiset('solver','ipopt','display','off');

% Initialization of return values and number of iterations for weighted sum
% approach to solve the subproblem
iterations = size(weights,2);
sol_x = zeros(n+m,iterations);
sol_f = zeros(p,iterations);

for i=1:iterations
    [sol_cont] = fmincon(@(x) weights(:,i)'*f([x;x_int]),x0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
%     opti_object = opti('fun',@(x) weights(:,i)'*f([x;x_int]),'ineq',Aineq,bineq,'eq',Aeq,beq,'nlmix',nonlcon,nonlcon_rhs,-ones(q,1),'bounds',lb,ub,'options',options);
%     sol_cont = solve(opti_object,x0);
    sol_f(:,i) = f([sol_cont;x_int]);
    sol_x(:,i) = [sol_cont;x_int];
end
end