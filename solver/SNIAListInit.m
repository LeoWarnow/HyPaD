function [int_lb,width,int_num,current_int]=SNIAListInit(n,m,lb,ub,SEED)
%SNIAListInit computes initial values for SNIAList

% Initialization
int_lb = lb(n+1:n+m);
int_ub = ub(n+1:n+m);
width = int_ub-int_lb;
int_num = prod(width+1);
% if SEED == 0
%     int_list = 0:int_num-1;
% else
%     rng(SEED, 'twister');
%     int_list = randperm(int_num)-1;
% end
current_int = 0;
end