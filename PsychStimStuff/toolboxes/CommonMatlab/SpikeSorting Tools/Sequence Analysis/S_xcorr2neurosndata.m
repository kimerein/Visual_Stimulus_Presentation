unitX = 9;
unitY = 15;
range =30;
binsize = 1e-3;
binperdata = .001;
figheader = 'Xcorr';

unit_spiketimesX = spikes{unitX};
unit_spiketimesY = spikes{unitY};
% Bin Data
binned_spikesX = bin2(unit_spiketimesX,binperdata);
binned_spikesY = bin2(unit_spiketimesY,binperdata);

cross_corr = xcorrmanual(binned_spikesX,binned_spikesY,range);

figure(3);
% a = log(cross_corr);
% b = find(a == -Inf); a(b) =0;
plot([1:size(cross_corr,2)]-(range+1),cross_corr/max(cross_corr),'-g ');% X axis time in Sec
axis tight
xlabel('ms');
%         temp = sprintf('%s%d-%d (%d,%d) Bin:%1.2gms',figheader,k,j,unit_Nspikes(k),unit_Nspikes(j),binsize*1000);
temp = sprintf('%s%d-%d Spikes:%d,%d, Bin:%1.2gms',figheader,unitX,unitY,1,1,binsize*1000);
title(temp);