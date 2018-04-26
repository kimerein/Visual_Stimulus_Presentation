function [params] = rmSkparams(skipped,params)
% function [params] = rmSkparams(skipped,params)
% 
% removed params that are skipped
% e.g. USAGE works for removing skipped by getWOI from findFracPeak parameters
%

notskipped = find(~skipped);

params(1).ind_insw = params(1).ind_insw(notskipped); %% MEM make into int32
params(1).risetime = params(1).risetime(notskipped) ;
params(1).eventsw = params(1).eventsw(notskipped); %% MEM make into int16
params(1).max2min = params(1).max2min(notskipped);
params(1).ind_fracpeak = params(1).ind_fracpeak (notskipped);
params(1).max = params(1).max (notskipped);
%  params(1).skipped = params(1).skipped (notskipped);