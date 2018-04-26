function [xEdge] = quartile(a,x)
%% function quartile(a,x)
% INPUT: histogram data
%                 frequency 
%                 value of bin
% OUTPUT: xEdge = bottom most bin value for each quadrent
%          i.e. xEdge(1) = bottom of 4th quadrant
%               xEdge(2) = bottom of 3rd quadrant etc..
%           length(xEdge) = 3 because first quadrant starts with the lowest
%           x value

l = 0.25; %% quadrant (not some other fraction of distribution
s = 0 ;
suma = sum(a);


%% take top 25perc
i = length(a);  %% get top l-fraction of distribution
while ( s < l*suma& i >0)
    s  = s+ a(i);
    i = i-1;
end
j=1;
xEdge(j)= x(i+1);
%% get 2nd largest l-fraction
s = 0;
while ( s < l*suma& i >0)
    s  = s+ a(i);
    i = i-1;
end
j=2;
xEdge(j)= x(i+1);
%% get 3rd largest l-fraction
s = 0;
while ( s < l*suma& i >0)
    s  = s+ a(i);
    i = i-1;
end
j=3;
xEdge(j)= x(i+1);

