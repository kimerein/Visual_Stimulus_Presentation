function xcorrmk(X,Y,allspikeTimes,ID,Ts)
% function xcorrmk(X,Y,allspikeTimes,label,Ts);
%
% DESC: crosscorrelated spikes with labels in X with those with labels in Y.
% used to cross-correlated isolated units or klusters
%
% INPUT:
% X - vector of labels (e.g. [1 2 4 5..]
% Y - vector of labels
% allspikeTimes - time of spikes in experiment (samples)
% ID - label for each of the spikes in allspikeTimes, size(ID)==size(allspikeTimes)
% Ts - sampling interval (sec)
%
% OUTPUT: length(X) figures each with length(Y) subplots
%         can also save to disk
%
% NOTE: xcorr is done manually and range is hardcoded
%                                    binsize is hardcoded
% BA 012307

srebin = 5;
range =500/srebin;

binsize = srebin*1e-3;
binperdata = Ts/1000;
figheader = 'KlustersXcorr';

NX = size(X,2);
NY = size(Y,2);
NUnits = NY;
NS = min(5,NUnits);
COL = ceil(NUnits/NS);
NN = NX;
if X == Y
    NN =     ceil(NX/2);
end
for k = 1:NN % Number of figures (1 for each cell)
    figure(k+700);
    close(k+700);
    figure(k+700);
    hold off
    set(gcf,'Position',[700 0 700 1120])
    
    unitX = X(k);
    unit_spiketimesX = allspikeTimes(find(ID == unitX));
    binned_spikesX = rebin(bin2(unit_spiketimesX,binperdata)',srebin);
   %% rebin normalization so that area is remains correct 
    for j = 1:NY %% Number of subplots per figure (1 for each cell)
    unitY =Y(j);
    unit_spiketimesY = allspikeTimes(find(ID == unitY));
        subplot(NS,COL,j)
       % Bin Data
        binned_spikesY = rebin(bin2(unit_spiketimesY,binperdata)',srebin);
        %         plot((1:size(binned_spikes,2))*(binsize),binned_spikes(1,:))  % X axis time in Se

        cross_corr = xcorrmanual(binned_spikesX,binned_spikesY,range);
        plot(([1:size(cross_corr,2)]-(range+1))*binsize*1000,cross_corr,'-b');% X axis time in Sec
%         plot(cross_corr,'-b');% X axis time in Sec
        axis tight
        xlabel('ms');
%         temp = sprintf('%s%d-%d (%d,%d) Bin:%1.2gms',figheader,k,j,unit_Nspikes(k),unit_Nspikes(j),binsize*1000);
        temp = sprintf('%s%d-%d Spikes:%d,%d*, Bin:%1.2gms',figheader,unitX,unitY,length(unit_spiketimesX),length(unit_spiketimesY),binsize*1000);
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