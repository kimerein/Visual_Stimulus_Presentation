function [binCenterx valuey weight] = binDataW(datax,datay,binsx)
% function [binCenterx valuey weight] = binDataW(datax,datay,binsx)
%
% INPUT: 
%        datax, datay: vectors of the same length.
%        binsx bin edges (see binIndex)
% OUTPUT: 
%        binCenterx vector containing the centers of datax bins length(binCenterx) = length(binsx) -1
%        valuey average value of datay in each bin, same length as
%        binCenterx
%        weight number of data points in each bin, same length as
%        binCenterx
%
% BA070608
[a binCenterx]= binIndex(datax,binsx);
for j=1:length(binCenterx)
    temp = a==j;
    weight(j) = sum(temp); % create a "weight" value that proportional to N,
    valuey(j) =  mean(datay(temp));
end