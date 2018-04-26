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


figure(1)       % rows X-Y (Yaxis), columns X-Z
imagesc(seq)
stemp = sprintf('%d-%d|%d-%d \n 0 = 31, %1.1f bins',unitX, unitZ, unitX,unitY,binsize*1000);
title(stemp);
stemp = sprintf('%d-%d',unitX,unitY);
ylabel(stemp);

% Project subset of data onto X-axi
test=   sum(seq(28:30,:),1);
figure(4)
plot([1:(range*2+1)]-(range+1),test,'k');


%Significance.
% randomUnit will randomly decided if there is an event in each time
% window (given N events in experiment of length T.
% then cross correlate
N = 177+35+20;%events
T= 1200000 ;% ms
nrepeats = 100;
randxcorr = zeros(nrepeats,2*range+1,'int32');
for i = 1:nrepeats;
    randomUnit = int32( rand(1,T) < N/T); % note in 1ms binned format
    randxcorr(i,:) = xcorrmanual(binned_spikesZ,randomUnit,range);
end
a= mean(mean(randxcorr,1));
b = mean(std(single(randxcorr),1));

figure(2); % plot 2Ds
plot(a+b*2);
axis tight;
stemp = sprintf('%d-RandUnit (95%% confidence)',unitX);
title(stemp);

%% calculate again assuming 9 random
   