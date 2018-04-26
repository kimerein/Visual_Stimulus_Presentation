function [index] = findvalue(data,value,mode)
% function [index] = findvalue(data,value,mode)
% returns index in data that contains number closes to value.
% if mode = -1 returns the first index found
%           =  1 returns the last index found
if nargin <3; mode = 0; end
switch(mode)
    case -1
        index = find(abs(data-value) == min(abs(data-value)),1,'first');
    case 1
        index = find(abs(data-value) == min(abs(data-value)),1,'last');
    otherwise
        index = find(abs(data-value) == min(abs(data-value)));

end

% intraindex = interp1(1:length(data),data,value);
% intrapolate
