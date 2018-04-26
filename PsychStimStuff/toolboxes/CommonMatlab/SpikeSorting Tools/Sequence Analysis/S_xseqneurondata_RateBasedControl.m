unitX = 7; unitY = 15; unitZ = 9;
range = 30;
binsize = 1e-3;
binperdata = .001;
      unit_spiketimesX = spikes{unitX};
    binned_spikesX = bin2(unit_spiketimesX,binperdata);
      unit_spiketimesY = spikes{unitY};
    binned_spikesY = bin2(unit_spiketimesY,binperdata);
      unit_spiketimesZ = spikes{unitZ};
    binned_spikesZ = bin2(unit_spiketimesZ,binperdata);
    seq = seqmanual(binned_spikesX,binned_spikesY,binned_spikesZ,range);
    figure(1)
    imagesc(seq)
    stemp = sprintf('%d-%d|%d-%d \n 0 = 31, %1.1f bins',unitX, unitZ, unitX,unitY,binsize*1000);
     title(stemp);
     stemp = sprintf('%d-%d',unitX,unitY);
   ylabel(stemp); 
    % rows X-Y (Yaxis), columns X-Z
    
    
% Project subset of data onto X-axi
test= sum(seq(28:30,:),1);
figure(4)
plot([1:(range*2+1)]-(range+1),test,'k');


%Significance.
% randomUnit will randomly decided if there is an event in each time
% window (given N events in experiment of length T.
% then cross correlate
% N = 177+35+20;%events
Ratewindow = 100;%ms
T= 1200000 ;% ms

%%% 
binned_spikesZ2 = bin2(unit_spiketimesZ,binperdata*Ratewindow); %% get rate in 100ms windows
binned_spikesZ2 = [binned_spikesZ2 zeros(1,T/Ratewindow - size(binned_spikesZ2,2))]; %% might not be spikes at the end of experiments so there will be missing bins
[cross_corr inX] = xcorrmanual(binned_spikesX,binned_spikesY,range);
binned_spikesXY([inX{32}' inX{33}' inX{34}']) = 1;

nrepeats = 50;
randUnit = zeros(1,T);
randxcorr = zeros(nrepeats,2*range+1,'int32');
for i = 1:nrepeats;
    
    for j = 1:T/Ratewindow
        randomUnit(1,(j-1)*Ratewindow+1:j*Ratewindow) = int32( rand(1,size((j-1)*Ratewindow+1:j*Ratewindow,2)) < double(binned_spikesZ2(j))/Ratewindow); % note in 1ms binned format
    end
    randxcorr(i,:) = xcorrmanual(binned_spikesXY,randomUnit,range);
end
a= mean(mean(randxcorr,1));
b = mean(std(single(randxcorr),1));

a= mean(randxcorr,1);
b = std(single(randxcorr),1);

figure(4); % plot 2Ds
plot([1:(range*2+1)]-(range+1),a+b*2);
axis tight;
stemp = sprintf('%d-RandUnit (95%% confidence)',unitX);
title(stemp);

%% calculate again assuming 9 random


[cross_corr inX] = xcorrmanual(binned_spikesX,binned_spikesY,range);