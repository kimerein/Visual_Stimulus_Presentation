

savefileheader1 = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Data Analysis\Figures\';
if strfind(savefilename,'\')
    filename = savefilename(max(strfind(savefilename,'\'))+1:end)
else
    filename = savefilename;
end
savepicfilename = strcat(savefileheader1,filename);


figID = 100;
%PLOT Spiketimes (concatanated sweeps)
%Objective 001
close(figID+10)
figure(figID+10)
binsize = 1000; %(ms)
binned = bin(spiketime(:,1),(N_samples/N_chn * N_sweeps)/ Ts/(binsize*1e-3))/(binsize*1e-3);
tmean = mean(binned);
tstdev = std(binned);
plot((1:size(binned,2))*(binsize*1e-3),binned(1,:))
stemp = sprintf('%s \nBinned (%3.1fms) # vs Time \nm = %2.1f std = %2.1f thres: %3.1f',processedFiles,binsize, tmean, tstdev, threshold);
xlabel('Time (s)');
ylabel('# Spikes');
title(stemp,'Interpreter','none');

%LINES for Drugs
if false
    % data in SWEEPS
    data = [50 83 143 282]
    temp = N_samples/Ts*data;
    for(i=1:size(data,2))
        line([temp(i) size(binned,2)*(binsize*1e-3)],[max(binned(1,:))*(1-.04*i) max(binned(1,:))*(1-.04*i)],'LineWidth',5,'Color','k');
    end
end

% MEAN and std from over DRug intrevals
if false
    figure(figID +100)
    % data in sec
    data = [100 400; 630 747;930 1200;1380 2500; 2720 2907]
    for(i=1:size(data,1))
        tmean(i) = mean(binned(1,data(i,1):data(i,2)));
        tstd(i) = std(binned(1,data(i,1):data(i,2)));
    end
    % TOP is STD
    bar([tmean;tstd]','stack')
    stemp = sprintf('%s \nBinned (%3.1fms)\nMean (BLUE) Std (RED)',processedFiles,binsize);
    title(stemp,'Interpreter','none');
% SAVE To TIFF
figure(figID +100)
temp = sprintf('%sP02.tif',savepicfilename(1:end-4));
print('-dtiffn',temp)

temp = sprintf('%s_pharmdata',savefilename(1:end-4));
  save(temp,...
        'tmean','tstd','-mat');

end


figure(figID+10)
% SAVE To TIFF
temp = sprintf('%sP01.tif',savepicfilename(1:end-4));
print('-dtiffn',temp)

