function  plotsgram(x,dt,fid,prange)
% function   plotsgram(x,dt,fid,prange)
% INPUT:
%     x - vector of data
%     dt - 1/ sampling rate
%     fid - figure handle for plot
%     prange - (optional) [lower upper] limits in Hz to plot spectrogram
%          default [0 60]
%  OUTPUT:
%
%  DESC:   plots spectrogram  in 30s bins of time 
% [S,F,T,P]  = spectrogram(double(x),1/dt*tbin,[],[],1/dt,'yaxis');
%
% BA030806
if nargin <4
prange = [0 length(x)*dt-1] ;% Hz
end
tbin = prange/10; %(s) 


figure(fid)
set(gca,'Color','none','XGrid', 'on')
 spectrogram(double(x),1/dt*tbin,[],[],1/dt,'yaxis');
title('Spectrogram')
ylim(prange)
caxis([(max(caxis)-(max(caxis)-min(caxis))/5) max(caxis)])