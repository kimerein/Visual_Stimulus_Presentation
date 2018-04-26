
b_ONEPLOT = 1;
figID = 100;
%PLOT Sweep  Threshold & Patched Cell
if b_ONEPLOT
    figure(figID)
    subplot(2,2,2)
else
    figure(figID)
    subplot(1,2,1)
end
temp = [1:10];% [1:size(data,2)];
plot(time_data,data(:,temp),'b');
hold all;
stemp = sprintf('%s\n SwFreq(Hz)=%3.0f \nThres= %2.1f Ts= %5.0fkHz',processedFiles,1/(N_samples/N_chn/Ts),threshold,Ts/1000 );
% ylim([min(data(5000:end-50,1))*1.5 max(data(5000:end-50,1))*1.5]);
% ylim([-150 50]);
xlim([0 max(time_data)]);
title(stemp,'Interpreter','none');
xlabel('Time(s)')
ylabel('Voltage')
% threshold
line([0 time_data(end,1)], [threshold threshold],'Color','k'); %(x1 x2), (y1 y2)
if b_ONEPLOT
    subplot(2,2,4)
else
    subplot(1,2,2)
end

% %%PLOT SPIKE WAVEFORMS
%TODO turn to histogram
%   plot criteraila
% figure;
for i=1:size(extracted_spikes,2)
    plot([1:size(extracted_spikes,1)]*1000/Ts,extracted_spikes(:,i));
    hold on
end
stemp = sprintf('%s \nSpike # %d\nthres: %2.2f max:%2.1f maxwidth(ms):%1.1f',processedFiles,Nextactedspikes,threshold,maxpeakamp,maxpeakwidth*1000);
title(stemp,'Interpreter','none');
xlim([0 size(extracted_spikes,1)*1000/Ts]);
% ylim([-150 50]);
ylabel('Voltage');
xlabel('Time (ms)');
%
%%%%%%%%%%%%%%%%%%%%%%55555555



% Objective 003
% Patched Cell
if b_ONEPLOT
    figure(figID)
    subplot(2,2,3)
else
    figure(figID+10)
end
plot(time_data, dataIC(:,2),'r');

%Find INTERVAL time of Depolarization  %Secondary IC channel
% temp = dataIC2(500:end,1)+ abs(mean(dataIC2(500:550,1))) > .5;  % 500pA; %skip beginning because of garbage (unknown why) % offset to zero
toffset = round(size(dataIC2,1)/4);  %injected current only can occur in 2nd half of trace
temp = dataIC2(toffset:end,1);
temp1 = diff(temp);
tmax = find(temp1==max(temp1)); tmin = find(temp1==min(temp1));
if isempty(tmax) || isempty(tmin)
    Interval = [0 0];
% elseif (tmin(1)-tmax(1)) ~= sum(temp(tmax:tmin+1)>(temp(tmax+1)/2)) %current is continuously high between tmax and tmin
    Interval = [0 0];
else
    Interval = [tmax tmin]+toffset;
end
%     Interval = [0 0];


if b_ONEPLOT
    figure(figID)
    subplot(2,2,1)
else
    figure(figID+4)
end
%  plot(1:size(spiketimeSw,1),spiketimeSw(:,3)/Ts) % PLOT ISI
if true % plot 
    %CALC & PLOT Mean #of spikes in sumspiketimeSw
    % sumspiketimeSw collapses all sweeps into one sweep for an average
    meanspiketime = cumsum(sumspiketimeSw)/N_sweepsExtracted;
    % meanspiketime = diff(meanspiketime);

    plot([1:size(meanspiketime,2)]/Ts,meanspiketime,'LineWidth',1)
    stemp = sprintf('%s\n Mean Sweep (n=%d)\n Thres= %2.1f',processedFiles,N_sweepsExtracted,threshold);
    title(stemp,'Interpreter','none');
    ylabel('Mean Spikes');
    xlabel('Time (s)');
    temp = max(meanspiketime)/(size(meanspiketime,2)/Ts);
    stemp = sprintf(' %2.1f  Total Spikes/s',temp);
    text(3* size(meanspiketime,2)/Ts/6,.9*max(meanspiketime), stemp);
    xlim([0 size(meanspiketime,2)/Ts]);
    %Interval LINES
    % line(Interval/Ts,[max(cumspiketime)-1 max(cumspiketime)-1],'Color','r','LineWidth',3); %HORIZONTAL
    line([Interval(1)/Ts Interval(1)/Ts],[min(meanspiketime)+1 max(meanspiketime)-1],'Color','r','LineWidth',0.5); %VERTICAL
    line([Interval(2)/Ts Interval(2)/Ts],[min(meanspiketime)+1 max(meanspiketime)-1],'Color','r','LineWidth',0.5); %VERTICAL
    hold on;
    % FIT data BEFORE STIM and PLOT
    if Interval(1) == 0
        temp = size(meanspiketime,2)/2;
    else
        temp = Interval(1);
    end
    temp = meanspiketime(1:temp-1); % extract before Interval
    temp1 = [1:size(temp,2)];
    p = polyfit(temp1,temp,1);
    temp1 = [1:size(meanspiketime,2)];
    temp = polyval(p,temp1);
    plot(temp1/Ts,temp,'Color','g');
    % UNITY line
    % line([0 size(meanspiketime,2)/Ts],[0 max(meanspiketime)],'Color','k','LineWidth',0.5);
end
if true %
    % PLOT MEAN spikes per second
    sigma = size(meanspiketime,2)/Ts * .005; % 1% of in units of samples
    x = [0:1/Ts:sigma*3];  %convolve with Gaussian for each sweep
    Gdist = 1/(sigma*sqrt(2*pi))*exp(-(x).^2/(2*sigma^2));
    %PLOT Instantaneous frequency
    diffmeanspiketime = [NaN diff(meanspiketime)]*Ts/10000; %units spikes/sec
    smoothdiffmean = conv(Gdist,diffmeanspiketime);
    plot([1:size(smoothdiffmean,2)]/Ts ,smoothdiffmean)
    hold on;
    xlim([0 size(meanspiketime,2)/Ts]);
    stemp = sprintf('%s\n Thres= %2.2f Sigma(ms) = %2.1f',processedFiles,threshold,sigma*1000);
    title(stemp,'Interpreter','none');
    line([Interval(1)/Ts Interval(1)/Ts],[min(smoothdiffmean)+1 max(smoothdiffmean)-1],'Color','r','LineWidth',0.5); %VERTICAL
    line([Interval(2)/Ts Interval(2)/Ts],[min(smoothdiffmean)+1 max(smoothdiffmean)-1],'Color','r','LineWidth',0.5); %VERTICAL

    % FIT data BEFORE STIM and PLOT
    if Interval(1) == 0
        temp = size(meanspiketime,2)/2;
    else
        temp = Interval(1);
    end
    temp1 = max(find(isnan(smoothdiffmean)==1)); % find first index that is not a NaN
    tmean = mean(smoothdiffmean(temp1+1:temp-1)); % extract before Interval
    tstdev = std(smoothdiffmean(temp1+1:temp-1));
    temp1 = [temp1+1:size(diffmeanspiketime,2)];
%     temp = polyval(p,temp1);
    line([0 size(diffmeanspiketime,2)/Ts],[tmean+tstdev tmean+tstdev],'LineStyle','--','Color','g','LineWidth',0.5); %VERTICAL
    line([0 size(diffmeanspiketime,2)/Ts],[tmean-tstdev tmean-tstdev],'LineStyle','--','Color','g','LineWidth',0.5); %VERTICAL
%     plot(temp1/Ts,temp,'Color','g');
end
%
temp = sprintf('%sP02.tif',savefilename(1:end-4));
print('-dtiffn',temp)

