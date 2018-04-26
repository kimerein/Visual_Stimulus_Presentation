binsize = 1000;
b = bin(a,ceil(max(a)/1000));
plot((1:size(b,2))*(binsize*1e-3),b(1,:))



%LINES for Drugs
if true
    % data in SWEEPS
    data = [4 15 24]
    temp = 30*data;
    for(i=1:size(data,2))
        line([temp(i) size(b,2)*(binsize*1e-3)],[max(b(1,:))*(0.5-.04*i) max(b(1,:))*(0.5-.04*i)],'LineWidth',5,'Color','k');
    end
end

% MEAN and std from over DRug intrevals
figID = 1;
if true
    figure(figID +100)
    % data in sec
    clear tmean; clear tstd
    data = [300 450; 550 700; 760 818]
    for(i=1:size(data,1))
        tmean(i) = mean(b(1,data(i,1):data(i,2)));
        tstd(i) = std(b(1,data(i,1):data(i,2)));
    end
    % TOP is STD
    bar([tmean;tstd]','stack')
    stemp = sprintf('%s \nBinned (%3.1fms)\nMean (BLUE) Std (RED)',processedFiles,binsize);
    title(stemp,'Interpreter','none');
% SAVE To TIFF
figure(figID +100)
temp = sprintf('%sP02.tif',savepicfilename(1:end-4));
print('-dtiffn',temp)
spiketimes = b;
temp = sprintf('%s_pharmdata.mat',savefilename(1:end-4));
  save(temp,...
       'spiketimes', 'tmean','tstd','-mat');

end


