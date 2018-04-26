%PLOT Spiketimes (concatanated sweeps)
%Objective 001
figure(figID+1)
binsize = 1000; %(ms)
binned = bin(spiketime(:,1),(N_samples/N_chn * N_sweeps)/ Ts/(binsize*1e-3))/(binsize*1e-3);
tmean = mean(binned);
tstdev = std(binned);
plot((1:size(binned,2))*(binsize*1e-3),binned(1,:))
stemp = sprintf('%s \nBinned (%3.1fms) # vs Time \nm = %2.1f std = %2.1f thres: %3.1f',processedFiles,binsize, tmean, tstdev, threshold);
xlabel('Time (s)');
ylabel('# Spikes');
title(stemp,'Interpreter','none');

% SAVE To TIFF
temp = sprintf('%sP01.tif',savefilename(1:end-4));
print('-dtiffn',temp)