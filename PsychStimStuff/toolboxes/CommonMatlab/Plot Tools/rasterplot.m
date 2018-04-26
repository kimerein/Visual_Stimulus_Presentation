function h = rasterplot(data,dt,swl,off,fid,len,scolor)
% function h = rasterplot(data,dt,swl,off,fid,len,scolor)
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
%           scolor - default 'k' character of 1x3 vector with color
% BA 020807

%%
if isempty(swl) || nargin < 3; swl = max(data)+1; end
if isempty(off) || nargin < 4; off = 10; end
if isempty(scolor) || nargin < 7; scolor = []; end
sw = floor(data./swl) ; % find which sweep a spike occurs in
data = data - sw*swl; % zero so that each spiketime is in reference to the correct episode (not beginning of the experiment)

if isempty(scolor) % default to use next color in colorOrder of 
scolor = [ 0 0 0];
end
% convert sw to trial (sweesp might be unevently spaced
% can't figure out how to do this without a loop
trial = nan(size(sw));
temp = unique(sw);
for ntrial = 1: length( temp)
    trial(sw==temp(ntrial))=ntrial;
end

%% plot
if isempty(len) || nargin<6; len = off/3; end
if isempty(fid) || nargin <5; fid=figure(); end
set(0,'CurrentFigure',fid);
% set(gcf, 'Renderer','openGL')]

if length(dt) ==1
       h = line(repmat(data',2,1)*dt*1000,repmat(trial',2,1)*off+repmat([-len/2;len/2],1,length(trial)),'linewidth',2,'color',scolor);
    
    %     for i = 1:length(data)
    %         line([data(i) data(i)]*dt*1000,[(sw(i)+1)*off-len/2 (sw(i)+1)*off+len/2],'linewidth',2,'color',scolor);
    %     end
else
    for i = 1:length(data)            
           h(i) =  line([data(i) data(i)]*dt(i)*1000,[(sw(i)+1)*off-len/2 (sw(i)+1)*off+len/2],'linewidth',2,'color',scolor);
    end
    set(fid,'Visible','on')
end
    %% old version
    
    % plot(data*dt,sw*off,'.b')
