function [out] = sem(data,dim)
%function [out] = sem(data,dim)
% computes standard error of the mean (SEM) = std/sqrt(n)
% dim specifies the dim along which to take the SEM.
% 1 means along a row (so the result is the std of the columns)
if nargin <2
dim = 1;
end
if dim ==1
    otherdim = 2;
else 
    otherdim = 1;
end
if isvector(data)
    if ~isrowvector(data)
        data = data';
    end
end
out = std(data,dim)/sqrt(size(data,otherdim));

