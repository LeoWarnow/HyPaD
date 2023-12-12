function [solutions,exitflag]=RSUP(n,m,p,q,f,g,Df,Dg,Aineq,bineq,Aeq,beq,lb,ub,T,llb,r)

clear model;

% Get sizes
sizeT = size(T,2);

% Initialize model
model.modelsense = 'min';
model.modelname = 'MasterProblem';

% Variables
model.vtype = [repmat('C',1,n),repmat('I',1,m),repmat('C',1,1),repmat('C',1,p)];
modelsize = n+m+1+p;

% Objective function
model.obj = [zeros(1,n+m),1,zeros(1,p)];

% Constraint function
sizeAineq = size(Aineq,1);
A = zeros(sizeAineq+sizeT*(p+q)+p,modelsize);
b = zeros(sizeAineq+sizeT*(p+q)+p,1);

% Original linear constraints
A(1:sizeAineq,:) = [Aineq, zeros(sizeAineq,1+p)];
b(1:sizeAineq) = bineq;

% Linearized constraints
for i=1:sizeT
    A(sizeAineq+((i-1)*p+1:i*p),:) = [Df(T(:,i)) zeros(p,1) -eye(p)];
    b(sizeAineq+((i-1)*p+1:i*p)) = -f(T(:,i))+Df(T(:,i))*T(:,i);
    A(sizeAineq+sizeT*p+((i-1)*q+1:i*q),:) = [Dg(T(:,i)) zeros(q,p+1)];
    b(sizeAineq+sizeT*p+((i-1)*q+1:i*q)) = -g(T(:,i))+Dg(T(:,i))*T(:,i);
end

% PS constraint
A(sizeAineq+sizeT*(p+q)+(1:p),:) = [zeros(p,n+m) -r eye(p)];
b(sizeAineq+sizeT*(p+q)+(1:p)) = llb;

% Equality constraints
if ~isempty(beq)
    sizeAeq = size(Aeq,1);
    A = [A;Aeq,zeros(sizeAeq,1+p)];
    b = [b;beq];
    sizeA = size(A,1);
    model.sense = [repmat('<',1,sizeA-sizeAeq),repmat('=',1,sizeAeq)];
else
    model.sense = '<';
end

model.A = sparse(A);
model.rhs = b;

% Upper and lower bounds
if ~isempty(ub)
    model.ub = Inf(1,modelsize);
    model.ub(1:n+m) = ub';
end
if ~isempty(lb)  
    model.lb = -Inf(1,modelsize);
    model.lb(n+m+1) = 0;
    model.lb(1:n+m) = lb';
end

% Warm start
% model.start = sol_start;

% Write model
gurobi_write(model, 'solver/RSUP.lp');

% Solve model and return results
clear params;
params.outputflag = 0;
result = gurobi(model, params);
if strcmp(result.status,'OPTIMAL')
    exitflag = 1;
    solutions = result.pool.xn;
elseif strcmp(result.status,'INFEASIBLE')
    exitflag = -1;
    solutions = -1;
else
    exitflag = -2;
    solutions = -2;
end
end