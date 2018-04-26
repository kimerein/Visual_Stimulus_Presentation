function [outind nomatch] = findSimEvent(eventx,eventy,window,range)
% function [outind nomatch] = findSimEvent(eventx,eventy,window,range)
% BA 100206
% script used in experiments where currents are recorded in two cells.
% 
% eventx and eventy are vectors containing event times
% window is the time in the same units as eventx and eventy that is to be
% considered "simultaneous"
%
% OUTPUT:
% outind is a 2 by n matrix.  where n is the length of x.
%       column 1 is the index of eventx and column 2 is the index of eventy that
%        occurs simultaneously
% nomatch - vector of indices of eventx that are unmatched to eventy
if nargin < 4,
    range = 1111111;
end
ii=0;
% out = [];
outind = [];
lastin = 1;
nomatch = [];
for i = 1: length(eventx)
    in = match(eventx(i),eventy(lastin:min(lastin+range,length(eventy))),window) ;

    if in~=-1

        in  = in +lastin -1;
%         ii = ii+1;
%         out(ii,:) = [eventx(i) eventy(in)];
        lastin = in;
        outind = [outind; [i in]];
    else
        in;
        nomatch = [nomatch i];
    end
    
end