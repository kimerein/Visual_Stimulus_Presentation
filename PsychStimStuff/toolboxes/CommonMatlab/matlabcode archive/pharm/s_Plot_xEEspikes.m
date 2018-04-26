
% readdirheader =  'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\';
% readdirheader =  '\\132.239.158.164\Scanziani Lab\Patch Data\';
readdirheader =  'e:\My Documents\Scanziani Lab\Patch Data\';
readfileheader = ['sr2005_10_21_'];
readfilenumber  = [2 3];
[filename] = createaxonfilename(readfileheader,readfilenumber(1));
pathfilename = strcat(readdirheader,filename);
[sfilename] = creatematlabfilename(readfileheader,readfilenumber(1));
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
        spiketime = spiketime +  spiketimeComb(size(spiketimeComb,1));
        N_sweepsExtracted = N_sweepsExtracted + N_sweepsExtractedComb(j-1);
    end
    % combine files
    processFilesComb = [processedFilesComb; processedFiles];
    N_sweepsComb = [N_sweepsComb;N_sweeps];
    N_sweepsExtractedComb = [N_sweepsExtractedComb; N_sweepsExtracted];
    N_samplesComb = [N_samplesComb;N_samples];
    spiketimeComb = [spiketimeComb;spiketime];
end
%BIN data using spike times (by sweep)
Ts = 100e-6; %[s/sample]
binwindow = 100e-3;%[s/bin]
samplesperBin = binwindow/Ts;
binData = zeros(N_sweepsExtractedComb(size(N_sweepsExtractedComb)),N_samplesComb/samplesperBin);
for(i = 1:size(spiketimeComb,1))
    sweepnum = 1+floor(spiketimeComb(i)/N_samples);
    binnum = 1+floor((spiketimeComb(i) - (sweepnum-1)*N_samplesComb)/samplesperBin);
    binData(sweepnum,binnum) =  binData(sweepnum,binnum) + 1;
end

%PLOT
figure;
% plot(binData(1:size(binData,1)*size(binData,2)))
plot(binData(1:104*size(binData,2)))
hold all;
% MARK FILES
fileBinnum = N_sweepsComb.* N_samplesComb/samplesperBin;

% MARK SWEEPS beginnings
clear sweeptime;
sweepnum = [33 93];
sweeptime = sweepnum*N_samples/samplesperBin;
sweeptime = [sweeptime; ones(1,size(sweeptime,2))]';
%     plot(sweeptime(:,1),sweeptime(:,2),'.r');

fileBinnum = (N_sweepsComb.* N_samplesComb/samplesperBin)';
fileBinnum = [fileBinnum; ones(1,size(fileBinnum,2))]';

plot(fileBinnum(:,1),fileBinnum(:,2),'.g');

