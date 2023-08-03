function [n,m,p,q,f,g,Df,Dg,Aineq,bineq,Aeq,beq,lb,ub,x0,is_convex,is_quadratic] = T6(~)
%T6 A non-quadratic test instance
%   This example was taken from:
%   Marianna De Santis, Gabriele Eichfelder, Julia Niebling, Stefan
%   Rocktäschel. Solving Multiobjective Mixed Integer Convex Optimization
%   Problems, SIAM Journal on Optimization (2020)

% Dimension of decision and criterion space
n = 2; % Continuous variables
m = 1; % Integer variables
p = 2; % Dimension criterion space
q = 1; % Number of constraints

% Problem type
is_convex = true;
is_quadratic = false;

% Objective function
f = @(x) [x(1)+x(3);x(2)+exp(-x(3))];
Df = @(x) [1,0,1;0,1,-exp(-x(3))];

% Linear constraints (Aineq*x <= bineq, Aeq*x = beq)
Aineq = [];
bineq = [];
Aeq = [];
beq = [];

% Lower and upper bounds (lb <= x <= ub)
lb = -2.*ones(n+m,1);
ub = 2.*ones(n+m,1);
% z = [-3;-1+exp(-2)];
% Z = [3;1+exp(2)];

% Start point x0
x0 = ceil((lb+ub)/2);

% Non-linear constraints (g(x) <= 0)
g = @(x) [x(1)^2+x(2)^2-1];
Dg = @(x) [2*x(1),2*x(2),0];
end