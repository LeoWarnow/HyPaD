function [n,m,p,q,f,g,Df,Dg,Aineq,bineq,Aeq,beq,lb,ub,x0] = P5(param)
%P5 A simple scalable convex example

% Dimension of decision and criterion space
n = param(2); % Continuous variables
m = param(2); % Integer variables
p = param(2); % Dimension criterion space
q = 2; % Number of constraints

% Objective function
f = @(x) [x(1:n)+x(n+1:n+m)];
Df = @(x) [eye(n),eye(m)];

% Linear constraints (Aineq*x <= bineq, Aeq*x = beq)
Aineq = [zeros(1,n),ones(1,m)];
bineq = [3];
Aeq = [];
beq = [];

% Lower and upper bounds (lb <= x <= ub)
lb = -3.*ones(n+m,1);
ub = 3.*ones(n+m,1);
% z = -6.*ones(p,1);
% Z = 6.*ones(p,1);

% Start point x0
x0 = ceil((lb+ub)/2);

% Non-linear constraints (g(x) <= 0)
g = @(x) [sum(x(1:n).^2)-1;sum(x(n+1:n+m).^2)-11];
Dg = @(x) [2.*x(1:n,1)',zeros(1,m);zeros(1,n),2.*x(n+1:n+m,1)'];
end