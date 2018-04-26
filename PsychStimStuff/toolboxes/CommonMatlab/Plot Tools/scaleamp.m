function sd = scaleamp(d,bmax)
% function sd = scaleamp(d,bmax)
% % scales d such that either the min or the max value in d is 1.
% INPUT: d - vector
%        bmax - optional (default<1>), if bmax = 0 then scaled to min
% BA0101208

if nargin <2
    bmax = 1;
end

sd = detrend(d);
if bmax
    sd = sd/max(sd);
else
    sd = sd/min(sd)*-1;
end