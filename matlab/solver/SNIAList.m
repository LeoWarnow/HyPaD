function [sol_x,flag,allSolved,current_int]=SNIAList(n,m,q,f,g,Aineq,bineq,Aeq,beq,lb,ub,intPNS,x_start,X,int_lb,width,int_num,current_int,a,r,OFFSET,GUESS)
%SNIAList Searches a new integer assignment by running through a list of
%all possible integer assignments

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
width_prod = prod(triu(ones(m),1).*(width'+1)+tril(ones(m)),2);
for i = current_int:int_num
%     number = int_list(i+1);
    number = i;
    x_int = int_lb+mod(floor(number./width_prod),width+1);
    current_int = i;
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