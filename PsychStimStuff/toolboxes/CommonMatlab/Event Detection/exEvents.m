function [outdata df]= exEvents(data,dt,fid,th,dfilter,bplot,params)
% function [outdata df]= exEvents(data,dt,fid,th,dfilter,bplot,params)
% INPUT: data = vector of data (remember to make sure the events are
% negative going i.e. IPSC should be flipped
%          dt   = 1/ sampling rate (s)
%          fid
%          th = (optional)string or number setting threshold for event detection
%                 e.g. th = 10; (same unit as data)
%               if th is a string it should be string containing a number.
%                 e.g. th = '1'; This the number is the number of standard
%                 deviations that should be used as the threshold.
%          dfilter =[HIGH_PASS_Hz LOW_PASS_Hz  boolean_60Hz_removal]
%                   e.g. [20 200 0]
%          params.dWOI
%                .frac
% OUTPUT:
%  outdata = structure w/ waveform, event statistics, threshold
%Desc: Function reads in VC/IC/LFP data  and extracts events from them
% alg:
% filter data
% thresholds events
% gets statistics of events
% gets waveform of events
% plots if bplot =1
% output data
%
% note filtering may not be optimal for all current types
%
% BA030907
tempfid = 3000;
ifid = tempfid;
if nargin < 4
    th = [];
end
if nargin < 6 || isempty(bplot)
    bplot = 1;
end
df = prepdata (data,bplot,dt,dfilter); %preprocess
if ischar(th)
    th = str2num(th)*mean(std(df));
end
[th ind_xth] = threshdata(df,dt,th); %threshold
if nargin <7
    WOId = [8 8]; % (samples) too look for first positive derivitive
    frac = 0.1; % fraction of peak to find
else
    WOId = params.WOId;
    frac = params.frac;
end

WOId = WOId.*round((1/dt/1000));

WOId = [50 30].*round((1/dt/1000)); % (samples) too look for first positive derivitive
frac = 0.5; % fraction of peak to find
outdata = findFracPeak(df,ind_xth,frac,WOId); % get parameters of event
% outdata = findFracPeakFilter(data, dt, dfilter, ind_xth,frac,WOId);
WOI =  [30 50].*round((1/dt/1000)); %window to extract
[event skipped] = getWOI(data,outdata.ind_fracpeak,WOI); %get waveform of event
outdata = rmSkparams(skipped,outdata);

outdata.th = th;
outdata.eventWV = event;
outdata.dfilter = dfilter;

bins = [1:5:1000]; % ms
if length(outdata.ind_fracpeak)>1
    [f x] = hist(diff(outdata.ind_fracpeak)*dt*1000,bins);
    try
        outdata.hIEI = [x; f];
    catch
        keyboard
    end
else
    outdata.hIEI = [bins; zeros(size(bins))];
end

%plots to check reasonableness of analysis
if bplot
    xtime = [1:size(outdata.eventWV,2)]*dt*1000;
    figure(ifid)
    clf
    tmp = min(size(outdata.eventWV,1),50)  % plot first 50 sweeps if there are more than 50
    plot(xtime,outdata.eventWV(1:tmp,:)-repmat(outdata.max(1:tmp)',1,size(outdata.eventWV,2)))
%     plot(xtime,outdata.eventWV)
    hold on
    plot(xtime,mean(outdata.eventWV),'-k','linewidth',2)
    axis tight
    title(num2str(size(outdata.eventWV,1)));
    
    ifid = infig(ifid);
    plot(x,f./sum(f))
    xlabel('ms')
    title('IEI')
%     xlim([max(1,x(find(f>1,1,'first'))) min(x(find(f>1,1,'last')+1),length(x))])
     xlim(xlimhist(x,f));
%REJECT based on IEI
end
bins = [5:10:2000]; % pA
[f x] = hist(outdata.max2min,bins);

if bplot
    ifid = infig(ifid);
    plot(x,f./sum(f))
    xlabel('pA/mV')
    title('pk2pk amp')  %% BA checked 10% of data pk2pk seems to work
    xlim(xlimhist(x,f));
    f = figs2subplots([tempfid:ifid],[ceil(sqrt(ifid-tempfid)) ceil(sqrt(ifid-tempfid))],[],fid);
    for i=tempfid:ifid
        close(i);
    end

end


