function [offset] = center(datain)
%%EXTRA
%% Center on Negative Peak 

if iscell(datain)
    for i =1 :size(datain,1)
        offset(i) = min(find(datain{i,1}(:,1,3) == min(datain{i,1}(:,1,3))));
    end
else
    for(i=1:size(datain,2))
        offset(i) = min(find(datain(:,i,3) == min(datain(:,i,3))));
    end
end