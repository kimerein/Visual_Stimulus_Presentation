function [outdata df]= exEventsF(data,dt,fid,th,dfilter,bplot,bth,WOI,params)
% function [outdata df]= exEventsF(data,dt,fid,th,dfilter,bplot,bth,WOI,params)
% INPUT: data = vector of data  (remember to make sure the events are
% negative going i.e. IPSC should be flipped)
% (assumed in mV)
%          dt   = 1/ sampling rate (s)
%          fid
%          th = (optional)string or number setting threshold for event detection
%                 e.g. th = 10; (same unit as data)
%               if th is a string it should be string containing a number.
%                 e.g. th = '1'; This the number is the number of standard
%                 deviations that should be used as the threshold.
%          dfilter =[HIGH_PASS_Hz LOW_PASS_Hz  boolean_60Hz_removal]
%                   e.g. [20 200 0]
%          bth = binary value 1 use specified threshold (default (0) sets a minimum to
%          threshold
%          WOI = [before after] window to extract in ms (default [30 30] 30 ms on each side of
%          event)
%          params.WOId % samples
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
% NOTE: similar to exEvents.m with the difference
%  1)that events are detected using the peak& trough detection and not thresholding
%  2) frindFracPeakFilter.m is used to find amplitude of unfiltered signal
%  (instead of findFracPeak.m which find the amplitude of the filtered
%  signal)
bdebug = 0;
tempfid = 3000;
ifid = tempfid;
if nargin < 4
    th = [];
end
if nargin < 6 || isempty(bplot)
    bplot = 1;
end
if nargin < 7 || isempty(bth)
    bth = 0;
end

df = prepdata (data,0,dt,dfilter); %preprocess
if ischar(th)
    th = str2num(th)*mean(std(df));
    if ~bth
        if abs(th) < 5; % assume usV
            th = -5;
        elseif abs(th)>30;
            th = -30;
        end
    end
end
if nargin <9
    WOId = [10 30]; %ms
    frac = 0.1; % fraction of peak to find
else
    WOId = params.WOId;
    frac = params.frac;
end

WOId = WOId.*round((1/dt/1000));% (samples) too look for first positive derivitive
if size(df,2)>size(df,1);
    df = df';
end
[ind_pk ind_tr] = detPeakTrough(df,th);
% warning('test BA'); th = th*.10;
% [th ind_xth] = threshdata(df,dt,th); %threshold
% ind_pk = ind_xth;
% 
% outdata = findFracPeakFilter(data, dt, dfilter, round((ind_pk+ind_tr)./2),frac,WOId);
if isrowvector(data)
    data = data';
end
outdata = findFracPeakFilter(data, dt, dfilter, ind_pk,frac,WOId);

if bdebug % plot detected events
    temp =data';
    hfid = plotdata(temp(1:end),dt);
    hold on;
    rasterplot(outdata.ind_fracpeak,dt,[],1e-4,hfid,[])
end
if nargin <8 || isempty(WOI)
    WOI =  [30 30].*round((1/dt/1000)); %window to extract
else
    WOI =  WOI*round((1/dt/1000)); %window to extract
end
% if isrowvector(data)
%     data = data';
% end
[event skipped] = getWOI(data,outdata.ind_fracpeak,WOI); %get waveform of event
outdata = rmSkparams(skipped,outdata);

outdata.th = th;
outdata.eventWV = event;
outdata.dfilter = dfilter;

bins = 30;%[1:5:1000]; % ms
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
    tmp = min(size(outdata.eventWV,1),50);  % plot first 50 sweeps if there are more than 50
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
% bins = [5:10:2000]; % pA
% [f x] = hist(outdata.max2min,bins);
[f x] = hist(outdata.max2min,20);

if bplot
    ifid = infig(ifid);
    plot(x,f./sum(f))
    xlabel('pA/mV')
    title('pk2pk amp')  %% BA checked 10% of data pk2pk seems to work
    %      xlim(xlimhist(x,f));

    f = figs2subplots([tempfid:ifid],[ceil(sqrt(ifid-tempfid)) ceil(sqrt(ifid-tempfid))],[],fid);
    for i=tempfid:ifid
        close(i);
    end

end


