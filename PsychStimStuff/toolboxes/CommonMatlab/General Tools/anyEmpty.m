function out = anyEmpty(data,dim)
% function out = anyEmpty(data,dim)
%
% function determines with each row/col of data contains an empty value
%
% INPUT: data: vector or matrix
%          dim (optional): search for rows (1, default) or columns(2)
%          containing empty values
%
% BA 062808 (UNTESTED)
if nargin <2; dim  = 1; end

if dim ==2
    data = data'; % change data to have the right shape for the selected dimension
end

out = zeros(size(data,1),1);
for i =1:size(data,1)
    if isempty( data(i,:))
        out(i) = 1;
    end
end

