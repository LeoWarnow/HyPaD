function [n,m,p,q,f,g,Df,Dg,Aineq,bineq,Aeq,beq,lb,ub,x0] = BL1(param)
%BL1 A test instance for BOMIP
%   This example was taken from:
%   Tyler Perini, Natashia Boland, Diego Pecin, Martin Savelsbergh. A
%   Criterion Space Method for Biobjective Mixed Integer Programming: The
%   Boxed Line Method, INFORMS Journal on Computing (2020)

% Dimension of decision and criterion space
n = 2; % Continuous variables
m = param+1; % Integer variables
p = 2; % Dimension criterion space
q = 1; % Number of constraints

% Randomized theta
rng(98693,'twister');
theta = [ones(1,param);zeros(1,param)];
theta_index = rand(1,param)>1/pi;
theta_rand_size = sum(theta_index);
theta(:,theta_index) = [1/4.*rand(1,theta_rand_size)+3/4;1/4.*rand(1,theta_rand_size)];

% Cones
k = 10;
d = k/(param+1)-0.5;
a = ones(1,param);
b = ones(1,param);
a(1) = -k+0.5*d;
b(1) = -a(1)-d;
for i=2:param
    a(i) = a(i-1)+2*d+1;
    b(i) = -a(i)-d;
end

% Objective function
f = @(x) [x(1);x(2)];
Df = @(x) [1,0,zeros(1,m);0,1,zeros(1,m);];

% Linear constraints (Aineq*x <= bineq, Aeq*x = beq)
Aineq = [-1,-1,zeros(1,param),2*k;-theta(1,:)',-(1-theta(1,:))',2*k.*eye(param),zeros(param,1);-theta(2,:)',-(1-theta(2,:))',2*k.*eye(param),zeros(param,1)];
bineq = [2*k;-theta(1,:)'.*a'-(1-theta(1,:))'.*b'+2*k*ones(param,1);-theta(2,:)'.*a'-(1-theta(2,:))'.*b'+2*k*ones(param,1)];
Aeq = [0,0,ones(1,m)];
beq = [1];

% Lower and upper bounds (lb <= x <= ub)
lb = [-k;-k;zeros(m,1)];
ub = [k;k;ones(m,1)];
% z = [];
% Z = [];

% Start point x0
x0 = ceil((lb+ub)/2);

% Non-linear constraints (g(x) <= 0)
g = @(x) 0;
Dg = @(x) zeros(1,n+m);
end