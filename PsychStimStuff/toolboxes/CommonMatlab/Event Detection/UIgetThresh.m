function ythres = UIgetThresh(data,dt)
% function ythres = UIgetThresh(data,dt)
%  note dt can be scalar (if sampling rate is the same throughout data)
%       dt vector of same length as data (if different sampling rate)
if length(dt)==1
    xtime = [1:length(data)]*dt;
else
    xtime = dt;
end
h = figure; %% PLOT to check threshold
clf
plot(xtime,data);
title('Click to SET threshold')
[x y] = ginput(1);
ythres = y;
close(h);