function [sol_x,B_new,flag,allSolved,current_int]=SNIA(n,m,p,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,intPNS,x_start,X,B,boxesLB,boxesUB,int_num,current_int,a,r)
%SNIA Searches a new integer assignment within predefined boxes

% Initialization
B_new = [];
if ~isempty(X)
    intArray = [intPNS{:,1},X(n+(1:m),:),B];
else
    intArray = [intPNS{:,1},B];
end
allSolved = false;

% Check if all integers are already solved
if ~(size(intArray,2)<int_num)
    allSolved = true;
    sol_x = lb;
    flag = 0;
    return;
end

% Find next integer assignment
boxes_num = size(boxesLB,2);
count = zeros(1,boxes_num);
for i=1:boxes_num
    count(i) = sum(all([boxesLB(:,i) <= intArray;intArray <= boxesUB(:,i)]));
end
[~,index] = min(count);
int_lb = boxesLB(:,index);
int_ub = boxesUB(:,index);
width = int_ub-int_lb;
width_prod = prod(triu(ones(m),1).*(width'+1)+tril(ones(m)),2);
for i = current_int(index):(prod(width+1)-1)
    x_int = int_lb+floor(i./width_prod);
    x_int(m) = int_lb(m)+mod(i,width(m)+1);
    current_int(index) = i;
    if m > 1
        if ~any(all(intArray==x_int))
            [sol_x,err]=FeasibilityCheck(n,q,g,Aineq,bineq,Aeq,beq,lb,ub,x_int,x_start);
            if err > 1e-6
                flag = 0;
            else
                flag = 1;
            end
            return;
        end
    else
        if ~any(intArray==x_int)
            [sol_x,err]=FeasibilityCheck(n,q,g,Aineq,bineq,Aeq,beq,lb,ub,x_int,x_start);
            if err > 1e-6
                flag = 0;
            else
                flag = 1;
            end
            return;
        end
    end
end
end