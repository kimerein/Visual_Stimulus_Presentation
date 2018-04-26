function [bindata] = rebin(data,binsperdata)
%[bindata] = bin(data,nbins)
% takes data vector of length L and puts ata in binsperdata into one bin
% bindata values are normalized by the number of databins that contributed
% to the bin

nbins =ceil(size(data,1)/binsperdata);
bindata = zeros(1,int32(nbins));
scal = zeros(1,int32(nbins));
for i=1:size(data,1) %%check right size
    bin = int32(ceil(i/binsperdata));
    bindata(bin) =  bindata(bin) + data(i);
    scal(bin) = scal(bin) + 1;
end
bindata = bindata./scal;