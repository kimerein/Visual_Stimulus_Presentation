function [bindata] = rebin(data,nbins)
%[bindata] = bin(data,nbins)
% takes data vector of length L and turns it into data vector of length
% NBins

bindata = zeros(1,int32(nbins));
 temp =size(data,1)/(nbins-1);
 
for i=1:size(data,1) %%check right size
    bin = i/nbins+1;
    bindata(bin) =  bindata(bin) + data(i);
end