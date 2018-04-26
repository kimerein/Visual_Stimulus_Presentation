function [evIM skippedIM] = getCAev(ev_ind,fdata,ftime,dtIM,fid)
%function [evIM skippedIM] = getCAev(ev_ind,fdata,ftime,dtIM,fid)
%  function retreives WOI around Ca_events and plots them 
%  INPUT 
%       ev_ind - vector of indices in fdata where events start
%       fdata - vector of Flou intensity
%       ftime - vector of frame times from begining of experiment (ms)
%              is used to interpolate between frames and create image data
%              at 1/dtIM
%       dtIM - (s) rate to interpolate IMdata at (usually 1/200Hz)
%       fid  - figure handle;
%
% BA 030407

tfINTER = interpIMdata(fdata,ftime,dtIM); % create interpolated flu signal
sel_xth_tfINTER = (ftime(ev_ind)/1000)/dtIM;% event time in realtime -> event time in index in tfINTER
WOI_evIM = round([.2 1]./dtIM); % sec
[evIM skippedIM] = getWOI(tfINTER',sel_xth_tfINTER,WOI_evIM);
figure(fid)
plot(evIM');
stemp = sprintf('N_{evCa}: %d',length(sel_xth_tfINTER));
title(stemp)