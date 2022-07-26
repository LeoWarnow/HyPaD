function [n,m,p,q,f,g,Df,Dg,Aineq,bineq,Aeq,beq,lb,ub,x0] = T7(~)
%T7 A quadratic test instance

% Dimension of decision and criterion space
n = 2; % Continuous variables
m = 2; % Integer variables
p = 2; % Dimension criterion space
q = 2; % Number of constraints


% Objective function
f = @(x) [x(1)+x(3);...
          x(2)+x(4)];
Df = @(x) [1,0,1,0;0,1,0,1];

% Linear constraints (Aineq*x <= bineq, Aeq*x = beq)
Aineq = [];
bineq = [];
Aeq = [];
beq = [];

% Lower and upper bounds (lb <= x <= ub)
lb = [-Inf;-Inf;-20;-20];
ub = [Inf;Inf;20;20];
% z = [];
% Z = [];

% Start point x0
x0 = [0;0;0;0];

% Non-linear constraints (g(x) <= 0)
g = @(x) [x(1)^2+x(2)^2-1;...
          (x(3)-2)^2+(x(4)-5)^2-10];
Dg = @(x) [2*x(1),2*x(2),0,0;0,0,2*(x(3)-2),2*(x(4)-5)];
end