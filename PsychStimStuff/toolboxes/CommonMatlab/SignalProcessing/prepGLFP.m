function df = prepGLFP (data,bplot,dt)
% function df = prepGLFP (data,bplot,dt)
%
% Preprocess Gamma LFP
% detrend, filter >1hz < 80Hz 
% BA030207
    sp = 'y';

if nargin < 3
    dt = 1;
end
if nargin < 2
    bplot = 1;
end
if bplot ==0
    sp = 'n';
end
% 
% params.Fs = 1/dt;
% params.fpass =  [1 100];
% params.tapers = [3 5];
% params.pad = 2;
% params.err = [1 0.0500];
% params.f0 = [];
% df =     rmlinesc(detrend(data),[],params,[],sp);
df = detrend(data);
df = filterdata(df,dt,20,1);
df = filterdata(df,dt,80,0);

if bplot
    figure(1)
    clf
plotdata(detrend(data,'constant'),dt,'trange',[0 3],'fid',1)
hold all
plotdata(df,dt,'trange',[0 3],'fid',1)
end