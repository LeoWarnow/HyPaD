function [sol_x,flag,allSolved,boxesLB,boxesUB]=SNIABoxesDynamic(n,m,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,intPNS,x_start,X,boxesLB,boxesUB,int_num,a,r,OFFSET,SELECTION_MODE,GUESS)
%SNIABoxesDynamic Searches a new integer assignment using a dynamic
%branching approach

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

% Find the box with least number of integer assignments
boxes_num = size(boxesLB,2);
count = zeros(1,boxes_num);
for i=1:boxes_num
    count(i) = sum(all([boxesLB(:,i) <= intArray;intArray <= boxesUB(:,i)]));
end
[count_min,index] = min(count);
branchLB = boxesLB(:,index);
branchUB = boxesUB(:,index);

% Compute the next boxes
if count_min > 0
    boxesLB(:,index) = [];
    boxesUB(:,index) = [];
    boxesLB(:,end+1) = branchLB;
    boxesUB(:,end+1) = branchUB;
    while true
        width = branchUB-branchLB;
        branchMids = branchLB+width./2;
        [~,k] = max(width); 
        countL = sum(intArray(k,:) < branchMids(k)+0.25);
        countU = sum(intArray(k,:) > branchMids(k)-0.25);
        if countL < countU
            branchUB(k) = floor(branchMids(k));
            boxesLB(k,end) = ceil(branchMids(k));
            boxesLB(:,end+1) = branchLB;
            boxesUB(:,end+1) = branchUB;
            if countL < 1
                break;
            end
        else
            branchLB(k) = ceil(branchMids(k));
            boxesUB(k,end) = floor(branchMids(k));
            boxesLB(:,end+1) = branchLB;
            boxesUB(:,end+1) = branchUB;
            if countU < 1
                break;
            end
        end
        intArray = intArray(:,all([branchLB-0.25 < intArray;branchUB+0.25 > intArray]));
    end
end

% Return next integer assignment
if SELECTION_MODE == 1
    x_int = branchLB;
elseif SELECTION_MODE == 2
    x_int = branchUB;
else
    x_int = round((branchLB+branchUB)./2);
end
[sol_x,err]=FeasibilityCheckPatch(n,q,g,Aineq,bineq,Aeq,beq,lb,ub,x_start,x_int);
if err < OFFSET
    flag = 1;
else
    flag = 0;
end
end