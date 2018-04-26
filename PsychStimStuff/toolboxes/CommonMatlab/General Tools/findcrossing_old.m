function [index] = findcrossing(data,value)
% function [index]= findcrossing(data,value)
% returns indices in data that are equal to or are closest to (the
% highside) of value
tmp = data-value;
index = find(diff(tmp>0));
