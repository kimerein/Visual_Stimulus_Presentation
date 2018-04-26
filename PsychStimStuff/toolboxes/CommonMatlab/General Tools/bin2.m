function [bindata] = bin2(data,binsperunit)
%[bindata] = bin2(data,binsperunit)
% This function histograms the data
% data is numbers describing an event.
% binsperunit is the number of bins for 1 unit of data.
% 
if isinteger(data)
    warning('function bindata does not accept Ints, data cast to single')
    data = single(data);
    binsperunit = single(binsperunit);
end
nbins =  floor((max(data))/binsperunit) +1 ;% First bin 0 >= x > binsperunit
                                           % Last bin takes last datapoint
bindata = zeros(1,nbins,'int32');
for i=1:length(data) %%check right size
    bin = floor(data(i)/binsperunit)+1;
    bindata(bin) =  bindata(bin) +1;
end