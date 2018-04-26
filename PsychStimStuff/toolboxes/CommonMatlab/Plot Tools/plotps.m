function [p1 f] = plotps(x,dt,fid,prange,param,pparam)
% function [p1 f] = plotps(x,dt,fid,prange,param,pparam)
% pparam.logplot = default 0
if nargin <3 || isempty(fid)
 fid=    figure;
else
    figure(fid);
end

if nargin <4 || isempty(prange)
    prange = [0 100];
end

if nargin <5; param = []; end

blogplot = 0;
if  exist('pparam')
    if isfield(pparam, 'logplot')
        blogplot = pparam.logplot;
    end
end

if min(size(x))==1 & ~isrowvector(x)
    x = x';
end
 
% WINDOW = length(x)/(1/dt)*10; % average in 1 sec windows
[p1 f param] = runps(x,dt,param);
if blogplot
    semilogy(f,p1)
    ylabel('log (units^2/Hz)')

else
    plot(f,p1)
    ylabel('(units^2/Hz)')

end
plotset(1);
xlim(prange)
title(['PS  win:' num2str(param.nwin*dt,'%1.2f') 's nfft:' num2str(param.nfft*dt,'%1.2f') 's']);
xlabel('Frequency (Hz)')
% resulution is half the sampling rate (1/dt/2) divided by NFFT.
% NFFT is the number of different frequency sinewaves used estimate the frequencis in the time series.
%
% Window function controls how you average over time.
% NFFT can be  greater than WINDOW (sharper)
% NFFT can be less than WINDOW (smoother)