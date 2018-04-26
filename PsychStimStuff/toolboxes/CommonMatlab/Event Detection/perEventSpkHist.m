function [eventHist skipped] = perEventSpkHist(eventtime, spiketime, swlength, WOI)
% function eventHist = perEventSpkHist(eventtime, spiketime, swlength, WOI)
% Performs a Peri-Event time histogram using event times as events.
%   INPUT:
%        eventtime -  times of events (in samples)
%        spiketime -  times of spikes (in samples)
%        swlength - length of sweep in samples
%        WOI - 1x2 array [before after] window of interest around eventtime
%        in samples
%   OUTPUT:
%       eventHist = spike times around each event
%           size(eventHist) = [length(eventtime) WOI(1)+WOI(2)+1]
%       skipped  = 1 if is skipped
bdebug = 0;
eventHist = zeros(length(eventtime),WOI(1)+WOI(2)+1,'int16');
xbin = [-WOI(1):WOI(2)];
skipped = zeros(length(eventtime),1,'int16');
ii=1;
for i=1:length(eventtime)
    % skip events that are within WOI of edge of sweep
    if floor((eventtime(i) - WOI(1))/swlength) == floor((eventtime(i) + WOI(2))/swlength)
        timediff = (spiketime-eventtime(i));
        inWOI = timediff(find((timediff >= -1 *WOI(1)) & (timediff <= WOI(2))));
        if bdebug
            if rem(i,50) == 0
                figure(10)
                plot(inWOI,'.b')
                title(num2str(find(abs(timediff)== min(abs(timediff)))))
                i
            end
        end
        [eventHist(ii,:) x] = hist(inWOI,xbin);
        
        ii = ii+1;
    else
        skipped(i) = 1;
    end
end