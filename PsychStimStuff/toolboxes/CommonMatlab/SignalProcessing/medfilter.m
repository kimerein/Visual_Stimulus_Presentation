function out = medfilter(data,window)
% out = medfilter(data,window)
% 
% INPUT : data - vector
%         window - 1 by 2 vector containing the number of steps (in
%         indicies) before and after a data point should be used when
%         calculating the median.
% OUPUT: out - median filtered data
% BA100607

ln = length(data);
out = zeros(size(data));
for i = 1:ln
    start = i-window(1);
    if  start<1;  start = 1;    end 

    stop = i+window(2);
    if  stop> ln;     stop = ln;    end
    out(i) = median(data(start:stop));
end