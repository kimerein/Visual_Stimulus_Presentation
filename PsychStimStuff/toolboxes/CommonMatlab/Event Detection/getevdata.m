function evdata = getevdata(data,dt,ind_xth)
% function evdata = getevdata(data,dt)
% Function extracts event waveforms and parameters
%
% INPUT: data - vector of data
%         dt  - 1/sampling rage
%         ind_xth  - index of threshold crossing
% OUPUT: evdata - structure containing event waveform, and parameters
%
% NOTE: waveform is extracted around ind_xth, not around ind_fracpeak
% BA 031007
df= prepGVC(data,1,dt);  %preprocess ( so that there is little baseline drift)
evdata.ind_xth = ind_xth;

WOId = [5 5].*round((1/dt/1000)); % (samples) too look for first positive derivitive
frac = 0.2; % fraction of peak to find
evdata = findFracPeak(df,ind_xth,frac,WOId); % get parameters of event

WOI =  [40 30].*round((1/dt/1000)); %window to extract
[event skipped] = getWOI(data,ind_xth(find(~evdata.skipped)),WOI); %get waveform of event
evdata = rmSkparams(skipped,evdata);

evdata.eventWV = event;
