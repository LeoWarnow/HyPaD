function [z,Z]=initBoundsNLP(n,m,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,x_int,x_start)
%initBoundsNLP Computes lower and upper bound for f(X_C,x_int)

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
options = optimoptions('fmincon','Display','none');

weights = eye(p);
z = zeros(p,1);
Z = zeros(p,1);
for i=1:p
    w = weights(:,i);
    [~,z(i)] = fmincon(@(x) w'*f([x;x_int]),x0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
end
Z = -Z;
end