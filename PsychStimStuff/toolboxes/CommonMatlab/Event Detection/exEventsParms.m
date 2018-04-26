function outdata = exEvents(df,dt,fid,th,bplot)
% function outdata = exEvents(data,dt,fid,th,bplot)
% INPUT: df = vector of data (remember to make sure the events are
% negative going i.e. IPSC should be flipped
%   NOTE: df  should be already preprocess/filtered
%          dt   = 1/ sampling rate (s)
%          fid
%          th = (optional)threshold for event detection
% OUTPUT:
%  outdata = structure w/ waveform, event statistics, threshold
%Desc: Function reads in VC/IC/LFP data  and extracts events from them
% alg:
% thresholds events
% gets statistics of events
% gets waveform of events
% plots if bplot =1
% output data
%
% note filtering may not be optimal for all current types
%
% BA030907
ifid = fid;
if nargin < 4
    th = [];
end
if nargin < 5 || isempty(bplot)
    bplot = 1;
end

[th ind_xth] = threshdata(df,dt,th); %threshold
WOId = [5 5].*round((1/dt/1000)); % (samples) too look for first positive derivitive
frac = 0.1; % fraction of peak to find
outdata = findFracPeak(df,ind_xth,frac,WOId); % get parameters of event
WOI =  [30 50].*round((1/dt/1000)); %window to extract
[event skipped] = getWOI(data,outdata.ind_fracpeak,WOI); %get waveform of event
outdata = rmSkparams(skipped,outdata);

outdata.th = th;
outdata.eventWV = event;

bins = [1:5:1000]; % ms
if length(outdata.ind_fracpeak)>1
    [f x] = hist(diff(outdata.ind_fracpeak)*dt*1000,bins);
    try
        outdata.IEI = [x; f];
    catch
        keyboard
    end
else
    outdata.IEI = [bins; zeros(size(bins))];
end

%plots to check reasonableness of analysis
if bplot
     ifid = fid +1;
    figure(ifid)
    clf
    stairs(x,f./sum(f))
    xlabel('ms')
    title('IEI')
    xlim([max(1,x(find(f>1,1,'first'))) min(x(find(f>1,1,'last')+1),length(x))])

    %REJECT based on IEI
end
bins = [5:10:2000]; % pA
[f x] = hist(outdata.max2min,bins);

if bplot
    ifid = ifid+1;
    figure(ifid)
    clf
    stairs(x,f./sum(f))
    xlabel('pA/mV')
    title('pk2pk amp')  %% BA checked 10% of data pk2pk seems to work
    xlim([max(1,x(find(f>1,1,'first'))) min(x(find(f>1,1,'last')+1),length(x))])
end
f = figs2subplots([fid:ifid],[ceil(sqrt(ifid-fid)) ceil(sqrt(ifid-fid))],[])
for i=fid+1:ifid
    close(i)
end


