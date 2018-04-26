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