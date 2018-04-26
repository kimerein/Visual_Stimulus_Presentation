selectUnits = [7 9 10 11 15];
NUnits = size(selectUnits,2);
range =30;
binsize = 1e-3;
binperdata = .001;
for k = 1:NUnits % Number of figures (1 for each cell)
    figure(k+500);
    hold off
    set(gcf,'Position',[700 0 700 1120])
    unitX = selectUnits(k);
    unit_spiketimesX = spikes{unitX};
    binned_spikesX = bin2(unit_spiketimesX,binperdata);
    for j = 1:NUnits %% Number of subplots per figure (1 for each cell)
        unitY = selectUnits(j);
        subplot(NUnits,1,j)
        unit_spiketimesY = spikes{unitY};
        % Bin Data
        binned_spikesY = bin2(unit_spiketimesY,binperdata);
        %         plot((1:size(binned_spikes,2))*(binsize),binned_spikes(1,:))  % X axis time in Se

        cross_corr = xcorrmanual(binned_spikesX,binned_spikesY,range);
%         plot(([1:size(cross_corr,2)]-(range+1))*binsize*1000,cross_corr,'-b');% X axis time in Sec
        plot(cross_corr,'-b');% X axis time in Sec
        axis tight
        xlabel('ms');
%         temp = sprintf('%s%d-%d (%d,%d) Bin:%1.2gms',figheader,k,j,unit_Nspikes(k),unit_Nspikes(j),binsize*1000);
        temp = sprintf('%s%d-%d Spikes:%d,%d, Bin:%1.2gms',figheader,unitX,unitY,1,1,binsize*1000);
        title(temp);
    end
    %% save data image to disk
    dirtemp = 'REPORT';
    figdesc = 'ManualXCorr';
    spath = [writedirheader dirtemp '\'];
    if (~isdir(spath))
        mkdir(spath);
    end
    temp = sprintf('%s%s_Unit_%d.png',spath,figdesc,unitX);
%     print('-dtiffn',temp);
    print('-dpng',temp);
end