function [boxesLB,boxesUB,int_num,current_int]=SNIABoxesInit(n,m,lb,ub,num_splits)
%SNIABoxesIniti Computes boxes to search for new integer assignments when
%using SNIABoxes

% Initialization
int_lb = lb(n+1:n+m);
int_ub = ub(n+1:n+m);
width = int_ub-int_lb;
int_num = prod(width+1);
boxesLB = int_lb;
boxesUB = int_ub;
for i=1:num_splits
    [max_width,k] = max(width);
    if max_width < 0.75
        break
    end
    offset = zeros(m,1);
    offset(k) = ceil(max_width/2-0.25);
    boxesLB = [boxesLB boxesLB+offset];
    boxesUB = [boxesUB-offset boxesUB];
    width(k) = floor(max_width/2+0.25);
end
current_int = zeros(1,size(boxesLB,2));
end