function df = prepGVC (data,bplot,dt,brm60)
% function df = prepGVC (data,bplot,dt,brm60)
%
% Preprocess Gamma VC reccording (detrend) >200
% so that there is little drift in baseline to allow effective thresholding
% of currents
% 
% BA030807
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

df = detrend(data);
if brm60
    params.Fs = 1/dt;
    params.fpass =  [1 100];
    params.tapers = [3 5];
    params.pad = 2;
    params.err = [1 0.0500];
    params.f0 = [];
    p = 10/length(data);
%     p = [];
    f0 = []; %60;
    df =     rmlinesc(detrend(data),f0,params,p,'y');
end
df = filterdata(df,dt,3,1);
df = filterdata(df,dt,200,0);


if bplot
    t = 5;% time in sec
    figure(1)
    clf
plotdata(detrend(data,'constant'),dt,'trange',[0 t],'fid',1)
hold all
plotdata(df,dt,'trange',[0 t],'fid',1)
% plotdata(df+std(df)*5,dt,'trange',[0 t],'fid',1)
end