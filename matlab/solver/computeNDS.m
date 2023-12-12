function [indexlist] = computeNDS(list)
%computeNDS Computes the nondominated points of a given set

size_list = size(list,2);
indexlist = true(1,size_list);

for i=1:size_list
    current_point = list(:,i);
    indexlist(i) = false;
    if ~any(all(current_point>=list(:,indexlist)))
        indexlist(i) = true;
    end
end
end