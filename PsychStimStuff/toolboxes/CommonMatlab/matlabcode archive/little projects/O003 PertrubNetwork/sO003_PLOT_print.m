%%%%%%%%%%%%%%%%%%%%%%
% Objective 003

figID = 100;

%RASTERPLOT %Spiketimes by Sweep
f = figure(figID+2);
subplot(1,2,1)
plot(spiketimeSw(:,2)/Ts,spiketimeSw(:,1),'.b');
stemp = sprintf('%s\n Thres= %2.2f',processedFiles,threshold);
title(stemp,'Interpreter','none');
ylabel('Sweep #');
xlabel('Time (s)');
%Interval
% line(Interval/Ts,[max(spiketimeSw(:,1))-1 max(spiketimeSw(:,1))-1],'Color','r','LineWidth',3); %HORIZONTAL
line([Interval(1)/Ts Interval(1)/Ts],[min(spiketimeSw(:,1))+1 max(spiketimeSw(:,1))-1],'Color','r','LineWidth',0.5); %VERTICAL
line([Interval(2)/Ts Interval(2)/Ts],[min(spiketimeSw(:,1))+1 max(spiketimeSw(:,1))-1],'Color','r','LineWidth',0.5); %VERTICAL
hold on;

% figure(figID+3)
subplot(1,2,2)
binsize = 30; %(ms)
binned = bin(spiketimeSw(:,2),N_samples/N_chn/Ts/(binsize*1e-3))/(binsize*1e-3);
tmean = mean(binned);
tstdev = std(binned);
plot((1:size(binned,2))*(binsize*1e-3),binned(1,:))
% line([Interval(1)/Ts Interval(1)/Ts],[min(binned)+1 max(binned)-1],'Color','r','LineWidth',0.5); %VERTICAL
% line([Interval(2)/Ts Interval(2)/Ts],[min(binned)+1 max(binned)-1],'Color','r','LineWidth',0.5); %VERTICAL
%% SOMETHING wrong with time scalling
line([0 size(binned,2)*binsize*1e-3],[tmean+tstdev tmean+tstdev],'LineStyle','--','Color','g','LineWidth',0.5); %HORizontal
line([0 size(binned,2)*binsize*1e-3],[tmean-tstdev tmean-tstdev],'LineStyle','--','Color','g','LineWidth',0.5); %HORizontal
line([0 size(binned,2)*binsize*1e-3],[tmean+2*tstdev tmean+2*tstdev],'LineStyle','--','Color','r','LineWidth',0.5); %HORizontal
line([0 size(binned,2)*binsize*1e-3],[tmean-2*tstdev tmean-2*tstdev],'LineStyle','--','Color','r','LineWidth',0.5); %HORizontal
line([0 size(binned,2)*binsize*1e-3],[tmean tmean],'LineStyle','--','Color','k','LineWidth',0.3); %HORizontal
stemp = sprintf('Binned (%2.1fms) \n%s\n Thres= %2.2f',binsize,filename,threshold);
title(stemp,'Interpreter','none');
ylabel('Spike #');
xlabel('Time (s)');
hold on
%% PLOT IC
if ~bpClamp
    plot((1:size(dataIC,1))/Ts,dataIC(:,1)/10+((tmean-tstdev) - max(dataIC(:,1)/10)),'Color','r')
    ylim([(tmean - 4*tstdev) (tmean + 4*tstdev)])
    xlim([0 max((1:size(dataIC,1))/Ts)])
end
%% SAVE TO EPS
% savefilename = sprintf('%seps',savefilename(1:end-3));
% print('-depsc2',savefilename)
% %SAVE TP PDF
% savefilename = sprintf('%spdf',savefilename(1:end-3));
% print('-dpdf',savefilename)
% SAVE To TIFF
% temp = sprintf('%sP01.tif',savefilename(1:end-4));
% print('-dtiffn',temp)