function xcorrk(X,Y,temp_spikeData,Ts);
% function xcorr2k(unitX,unitY,temp_spikeData,Ts);

range =200;
binsize = 1e-3;
binperdata = Ts/1000;
figheader = 'KlustersXcorr';

NX = size(X,2);
NY = size(Y,2);
NUnits = NY;
NS = min(5,NUnits);
COL = ceil(NUnits/NS);
for k = 1:NX % Number of figures (1 for each cell)
    figure(k+700);
    close(k+700);
    figure(k+700);
    hold off
    set(gcf,'Position',[700 0 700 1120])
    
    unitX = X(k);
    unit_spiketimesX = temp_spikeData{unitX}(:,1);
    binned_spikesX = bin2(unit_spiketimesX,binperdata);
    for j = 1:NY %% Number of subplots per figure (1 for each cell)
    unitY =Y(j);
    unit_spiketimesY = temp_spikeData{unitY}(:,1);
        subplot(NS,COL,j)
       % Bin Data
        binned_spikesY = bin2(unit_spiketimesY,binperdata);
        %         plot((1:size(binned_spikes,2))*(binsize),binned_spikes(1,:))  % X axis time in Se

        cross_corr = xcorrmanual(binned_spikesX,binned_spikesY,range);
        plot(([1:size(cross_corr,2)]-(range+1))*binsize*1000,cross_corr,'-b');% X axis time in Sec
%         plot(cross_corr,'-b');% X axis time in Sec
        axis tight
        xlabel('ms');
%         temp = sprintf('%s%d-%d (%d,%d) Bin:%1.2gms',figheader,k,j,unit_Nspikes(k),unit_Nspikes(j),binsize*1000);
        temp = sprintf('%s%d-%d Spikes:%d,%d, Bin:%1.2gms',figheader,unitX,unitY,1,1,binsize*1000);
        title(temp);
    end
%     %% save data image to disk
%     dirtemp = 'REPORT';
%     figdesc = 'ManualXCorr';
%     spath = [writedirheader dirtemp '\'];
%     if (~isdir(spath))
%         mkdir(spath);
%     end
%     temp = sprintf('%s%s_Unit_%d.png',spath,figdesc,unitX);
% %     print('-dtiffn',temp);
%     print('-dpng',temp);
end