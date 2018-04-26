function [offset] = centerIntra(datain)
%%EXTRA
%% Center on Maxslope

if iscell(datain)
    for i =1 :size(datain,1)
        offset(i) = find(diff(datain{i,1}(:,1,1)) == max(diff(datain{i,1}(:,1,1)))) +1;%add 1 since derivative reduces index
    end
else
    for(i=1:size(datain,2))
        offset(i) = find(diff(datain(:,i,1)) == max(diff(datain(:,i,1)))) +1 ;
    end
end