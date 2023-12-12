function [sol_x,flag,allSolved,current_int]=SNIABoxesFixed(n,m,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,intPNS,x_start,X,boxesLB,boxesUB,int_num,current_int,a,r,OFFSET,GUESS)
%SNIABoxes Searches a new integer assignment within predefined boxes

% Initialization
if ~isempty(X)
    intArray = [intPNS{:,1},X(n+(1:m),:)];
else
    intArray = [intPNS{:,1}];
end
allSolved = false;

% Check if all integers are already solved
if ~(size(intArray,2)<int_num)
    allSolved = true;
    sol_x = lb;
    flag = 0;
    return;
end

% Buest guess based on relaxation
if GUESS
    [x_relaxed,~,~,~] = PS(f,g,Aineq,bineq,Aeq,beq,lb,ub,x_start,a,r);
    sol_x_int = round(x_relaxed(n+1:n+m));
    [sol_x,err]=FeasibilityCheckPatch(n,q,g,Aineq,bineq,Aeq,beq,lb,ub,x_relaxed,sol_x_int);
    if err > 1e-6
        flag = 0;
    else
        flag = 1;
    end
    if m > 1
        if ~any(all(intArray==sol_x_int))
            return;
        end
    else
        if ~any(intArray==sol_x_int)
            return;
        end
    end
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
    x_int = int_lb+mod(floor(i./width_prod),width+1);
    current_int(index) = i;
    if m > 1
        if ~any(all(intArray==x_int))
            [sol_x,err]=FeasibilityCheckPatch(n,q,g,Aineq,bineq,Aeq,beq,lb,ub,x_start,x_int);
            if err < OFFSET
                flag = 1;
            else
                flag = 0;
            end
            return;
        end
    else
        if ~any(intArray==x_int)
            [sol_x,err]=FeasibilityCheckPatch(n,q,g,Aineq,bineq,Aeq,beq,lb,ub,x_start,x_int);
            if err < OFFSET
                flag = 1;
            else
                flag = 0;
            end
            return;
        end
    end
end
end