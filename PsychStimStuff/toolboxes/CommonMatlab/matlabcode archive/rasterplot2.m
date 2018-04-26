function rasterplot2(data,data2,dt,swl,off,fid,len)
% function rasterplot2(data,dt,swl,off,fid,len)
% creates rasterplot of spiketime data assuming a episode/sweep of
% swl(samples)
% 
% INPUT:
%        data - Vector of spike times ( time in samples)
%        dt - 1/sample rate (s)
%        swl - swlength (samples) <length of data>
%        off - offset along the y-axis  <10>
%        fid - figure handle <new figure>
%        len - <off/10>
%
% BA 020807
%%
if isempty(swl) || nargin < 3; swl = max(data)+1; end
if isempty(off) || nargin < 4; off = 10; end
sw = floor(data./swl) ; % find which sweep a spike occurs in
data = data - sw*swl; % zero so that each spiketime is in reference to the correct episode (not beginning of the experiment)

data2 = (data2 - min(data2))./max(data2) ; %normalize to be in 0-1 range
%% plot
if isempty(len) || nargin<6; len = off/3; end
if isempty(fid) || nargin <5; fid=figure(); end
figure(fid)
if length(dt) ==1
    for i = 1:length(data)
        line([data(i) data(i)]*dt*1000,[(sw(i)+1)*off-len*data2(i)/2 (sw(i)+1)*off+len*data2(i)/2],'linewidth',2,'color','k');
    end
else
    for i = 1:length(data)
        line([data(i) data(i)]*dt(i)*1000,[(sw(i)+1)*off-len*data2(i)/2 (sw(i)+1)*off+len*data2(i)/2],'linewidth',2,'color','k');
    end
end
%% old version

% plot(data*dt,sw*off,'.b')
