function outdata = exEvents(data,dt,fid,th)
% function outdata = exEvents(data,dt,fid)
% INPUT: data = vector of data (remember to make sure the events are
% negative going i.e. IPSC should be flipped
%          dt   = 1/ sampling rate (s)
%          fid
%          th = (optional)threshold for event detection
% OUTPUT:
%  outdata = structure w/ waveform, event statistics, threshold
%Desc: Function reads in VC/IC/LFP data  and extracts events from them
% alg:
% applies filtering
% thresholds events
% gets statistics of events
% gets waveform of events
% plots if bplot =1
% output data
%
% note filtering may not be optimal for all current types
%
% BA030907

if nargin < 4
    th = [];
end
bplot = 1;
%
df= prepGVC(data,1,dt);  %preprocess ( so that there is little baseline drift)

% clear data;

[th ind_xth] = threshdata(df,dt,th); %threshold
WOId = [5 5].*round((1/dt/1000)); % (samples) too look for first positive derivitive
frac = 0.2; % fraction of peak to find
CURparams = findFracPeak(df,ind_xth,frac,WOId); % get parameters of event

WOI =  [20 40].*round((1/dt/1000)); %
[event skipped] = getWOI(data,CURparams.ind_fracpeak,WOI); %get waveform of event
CURparams = rmSkparams(skipped,CURparams);

%plots to check reasonableness of analysis
if bplot
    xtime = [1:size(event,2)]*dt*1000;
    figure(fid)
    clf
%     plot(xtime,event)
    hold on
    plot(xtime,mean(event),'-k','linewidth',2)
    axis tight
    title(num2str(size(event,1)));
end
outdata.th = th;
outdata.eventWV = event;
outdata.ind_fracpeak = CURparams.ind_fracpeak;

bins = [1:100]; % ms
[f x] = hist(diff(CURparams.ind_fracpeak)*dt*1000,bins);
outdata.IEI = [x f];

if bplot
    figure(fid+1)
    clf
    stairs(x,f./sum(f))
    xlabel('ms')
    title('IEI')
    %REJECT based on IEI
end
bins = [5:10:2000]; % pA
[f x] = hist(CURparams.max2min,bins);
outdata.max2min = [x f];

if bplot
    figure(fid +2)
    clf
    stairs(x,f./sum(f))
    xlabel('pA')
    title('pk2pk amp')  %% BA checked 10% of data pk2pk seems to work
end





