function [n,m,p,q,f,g,Df,Dg,Aineq,bineq,Aeq,beq,lb,ub,x0] = P1(~)
%P1 A quadratic test instance
%   This example was taken from:
%   Kristo Mela, Juhani Koski and Risto Silvennoinen. Algorithm for
%   Generating the Pareto Optimal Set of Multiobjective Nonlinear
%   Mixed-Integer Optimization Problems, AIAA (2007)

% Dimension of decision and criterion space
n = 2; % Continuous variables
m = 8; % Integer variables
p = 2; % Dimension criterion space
q = 1; % Number of constraints


% Objective function
G = [1,-1,2,0,0,0,0,0,0,0;...
    -1,2,0,0,2,0,0,0,0,0;...
    0,0,3,0,2,0,0,0,0,0;...
    2,0,0,4,0,2,0,2,0,0;...
    0,0,0,0,5,2,0,0,0,0;...
    0,0,0,0,0,6,0,0,0,0;...
    0,0,0,0,0,0,7,0,0,0;...
    0,0,0,0,0,0,0,0,0,0;...
    0,0,0,0,0,2,0,0,0,0;...
    0,2,0,0,0,0,0,0,0,10];
c1 = [-1;-1;1;-10;0;1;-2;0;3;0];
c2 = [1;2;-1;1;5;-2;0;6;0;-3];
f = @(x) [0.5*x'*G*x+c1'*x;c2'*x];
Df = @(x) [(0.5*(G+G')*x+c1)';c2'];

% Linear constraints (Aineq*x <= bineq, Aeq*x = beq)
Aineq = [];
bineq = [];
Aeq = [];
beq = [];

% Lower and upper bounds (lb <= x <= ub)
lb = [-1;-1;0;0;0;0;0;0;0;0];
ub = ones(n+m,1);
% z = [-10;-10];
% Z = [25;15];

% Start point x0
x0 = ceil((lb+ub)/2);

% Non-linear constraints (g(x) <= 0)
g = @(x) 0;
Dg = @(x) zeros(1,n+m);
end