function [ind_xthres] = threshdiff(data, threshold)
% function [ind_xthres] = threshdiff(data, threshold)
% finds point of where the derivative of data is greater then threshold
%ind_xthres = find(diff(diff(data) > threshold)==1)+2;

ind_xthres = find(diff(diff(data) > threshold)==1)+2;