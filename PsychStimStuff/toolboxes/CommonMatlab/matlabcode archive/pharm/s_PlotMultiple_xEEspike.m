%SCRIPT to combined MAT files into continuous sweep
% READS MAT files with spiketime
% PLOTs cum spike number vs spike time
% PLOTs sweep number (entered by user)
clear all;

readdirheader =  'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\';
% readdirheader =  '\\132.239.158.164\Scanziani Lab\Patch Data\';
% readdirheader =  'e:\My Documents\Scanziani Lab\Patch Data\';
readfileheader = ['s2_2005_10_21_'];
savefileheader = ['slice1_2005_10_21.mat']
readfilenumber  = [2 3 4];
[filename] = createaxonfilename(readfileheader,readfilenumber(1));
pathfilename = strcat(readdirheader,filename);
[sfilename] = creatematlabfilename(savefileheader,readfilenumber(1));
savefilename = strcat(readdirheader,sfilename);

processedFilesComb = [];
N_sweepsComb = [];
N_sweepsExtractedComb = [];
N_samplesComb = [];
spiketimeComb = [];

%readfilenumbers
% plot in a row
binnumlastfile = 0;
%LOAD and COMBINE DATA
for (j = 1: size(readfilenumber,2))
    [filename] = creatematlabfilename(readfileheader,readfilenumber(j));
    filename =     strcat(readdirheader,filename);
    load(filename,...
        'processedFiles','spiketime','N_sweeps','N_sweepsExtracted','N_samples','Ts','subsample','-mat');
    if j >1
        N_sweepsExtracted = N_sweepsExtracted + N_sweepsExtractedComb(j-1);
        spiketime = spiketime +   N_sweepsExtractedComb(j-1)*N_samplesComb(j-1)/lastsubsample;
    end
    % combine files
    processFilesComb = [processedFilesComb; processedFiles];
    N_sweepsComb = [N_sweepsComb;N_sweeps];
    N_sweepsExtractedComb = [N_sweepsExtractedComb; N_sweepsExtracted];
    N_samplesComb = [N_samplesComb;N_samples];
    spiketimeComb = [spiketimeComb;spiketime];
    lastsubsample = subsample
end

save(sfilename,'spiketimeComb','-mat');
figure;
plot(spiketimeComb/Ts,[1:size(spiketimeComb,1)]);
xlabel('Time(s)');
ylabel('Cum Spikes');
stemp = sprintf('%s\nSpike # vs Spike Time',savefileheader)
title(stemp);
hold all;

% MARK SWEEPS beginnings
clear sweeptime;
sweepnum = [50 104+33 104+93 282]; %Slice1
sweeptime = sweepnum*90000/Ts;
sweeptime = [sweeptime; ones(1,size(sweeptime,2))]';
for i=1:size(sweepnum,2)
line([sweeptime(i) sweeptime(i)],[1+1 size(spiketimeComb,1)-1],'Color','r','LineWidth',0.5); %VERTICAL
end
% 
% fileBinnum = (N_sweepsComb.* N_samplesComb)';
% fileBinnum = [fileBinnum; ones(1,size(fileBinnum,2))]';

% plot(fileBinnum(:,1),fileBinnum(:,2),'.g');

