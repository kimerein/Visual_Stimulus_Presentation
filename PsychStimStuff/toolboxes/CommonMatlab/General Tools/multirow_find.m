function ind = multirow_find(data,findstuff)
% function ind = multirow_find(data,findstuff)
% finds the index of findstuff(i) in the i-th row of data
%
% BA 051508
ind = [];
j=1;
for i = 1:length(findstuff)
    temp = find(data(i,:)==findstuff(i));
    if ~isempty(temp)
    ind(j) = temp;
    j = j+1;
    end
end