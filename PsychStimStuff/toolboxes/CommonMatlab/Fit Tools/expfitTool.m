function s = expfitTool(data,dt,dwin)
% function fit_tau = expfitTool(data,dt,dwin)
% Fit data with exponential from over user selected interval
%                  s(1)+s(2)*exp(-tt/s(3));
% INPUT:
% data -  vector or matrix containing data to fit. If matrix then new piece of data in each row
% dt - sample interval
% dwin - (optional) window in data to use for fit (samples)
% OUTPUT:
%  s= exp fit units of dt
%  BA031809
%
figure(99);clf;
plotdata(data,dt,'fid',99,'color','.b'); axis tight
hold on;
if nargin <3; % get dwin
    a = ginput(2);
    dwin = round(a(:,1)/dt/1000);
end


tt = [dwin(1):dwin(2)]*dt;
s = exp2fit(tt,data(dwin(1):dwin(2)),1);
% plot fit  % from exp2fit example
fun = @(s,t) s(1)+s(2)*exp(-tt/s(3));
ff=fun(s,tt);
plot(tt*1000,ff,'r');
title(['Tau = ' num2str(s(3)*1e3)])

s = s*1e3; % dt is in sec but everything else is in ms