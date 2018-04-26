% Multielectrode waveform extraction
%   For each electrode threshold and extract spiketimes.
%   Collect waveform from all other electrodes (allow specifying how many
%   and how many samples)

%   Plot the raw data.
%   kmeans cluster
s_defineColors 

warning off
writedirheader = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Data Analysis\Expts\';
%%Data output Parameters
breport =0; % toggle report output
bsave = 0;  % toggle saving data e.g Xspikes,clusters..
b_intraspikes = 0;
bsortAmp =0;

% readdirheader =  'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\';
% [readfileheader readfilenumber]=  readaxonDialog(bSliceCell)
% processedFiles = createaxonfilename(readfileheader,readfilenumber(1))
% tic
% % output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\subset2006_03_20_0009.abf',-1,1);
% output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\test1.abf',-1,1,100,[2:8]);
output2 = [];
output1 = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_05_03_0006s.abf',-1,1,100,[2:8]);
lastINSW = output1.Nsweeps;%Last at inhibitory reversal sweep

output2 = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_05_03_0007s.abf',-1,1,100,[2:8]);


% output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_05_08_0003s.abf',-1,1);
% toc
%%% USER INPUT Data parameters
nIntracellsignals =1; % #ADC channels which were Intracellular (Code assumes they are contiguous from 1 to N)
indexoffset = 1  ; % #Offset to ADC channels cause data starts with time column
% neighbors = [9 10 11; 8 10 11; 8 9 11; 8 9 10] % ADC number electrodes that are used with other electrode in spikesorting (waveforms will be concatinated)
neighbors = [13 14;12 14;12 13]; % ADC number electrodes that are used with other electrode in spikesorting (waveforms will be concatinated)
% neighbors = [9 10 11 12;12 14;12 13]; % ADC number electrodes that are used with other electrode in spikesorting (waveforms will be concatinated)
EEgroup = [6 7 8];
% EEgroup = [2 3 4 5];
% neighbors = -1;

% ADD REPORT extracted
figure(1)
plot(output.data(:,2:output.Nchan:100));


%% ADD Preprocessing Filtering


%% CREATE DIR for data
writedirheader = [writedirheader output1.sfilename(1:end-4) '\' 'EE' strrep(num2str(EEgroup),'  ','_') '\'];
if ~isdir(writedirheader)
    mkdir(writedirheader);
end;

figheader = '';
%%% EXTRACTED from DATA
Ts = 1/output1.dt; % Sampling Frequency (Typically 20kHz)
N_chn = output1.Nchan; % Number of ADC Channels recorded from
N_samples = output1.Nsamples; % Number of ADC Channels recorded from
sourcefile = output1.sfilename;
sw = max(output1.xSweeps); % sweeps to set max time
sourcefile = output1.sfilename;

tempIND = min(size(output1.data,1),size(output2.data,1));
for i=1:nIntracellsignals
    IntraCellData = [output1.data(1:tempIND,1+i:output1.Nchan:end) output2.data(1:tempIND,1+i:output2.Nchan:end)];
end

%%% Spike Extraction Parameters
WOI = [400e-6 400e-6; 1100e-6 1100e-6];%Window of interest around Beg_spike -1000us + 2000us
std_threshold =4; %number of std for threshold
b_intra = 0; % trigger on negative going spike
maxpeakamp = 10000*std_threshold;  %in std
maxpeakwidth =  2e-3;


index_beg_spike = 0; temp1 = []; temp2 = []; temp2 = []; temp = []; fstd = []; fmean = [];
% for triggerEE=nIntracellsignals+1:N_chn  % all Extracellular signals , threshold
for i = 1: size(EEgroup,2)  % all Extracellular signals in 1 tetrode
    triggerEE = EEgroup(i);
    if ~isinteger(output1.data)
        fstd(i) =  mean(std([output1.data(1:tempIND,indexoffset+triggerEE:N_chn:5*N_chn) output2.data(1:tempIND,indexoffset+triggerEE:N_chn:5*N_chn)]));
    else
        fstd(i) =  mean(int16std([output1.data(1:tempIND,indexoffset+triggerEE:N_chn:5*N_chn) output2.data(1:tempIND,indexoffset+triggerEE:N_chn:5*N_chn)]));
    end

    fmean(i) =  mean(mean([output1.data(1:tempIND,indexoffset+triggerEE:N_chn:5*N_chn) output2.data(1:tempIND,indexoffset+triggerEE:N_chn:5*N_chn)]));

    %%% CHANGE code so spikes are in std units
    tic
    %% REmember threshold is mean - std*std_threshold
    temp = [temp; int32(FindSpikes2([output1.data(1:tempIND,indexoffset+triggerEE:N_chn:end) output2.data(1:tempIND,indexoffset+triggerEE:N_chn:end)], (fmean(i)-(std_threshold*fstd(i))) ,Ts,b_intra));];
    toc
    nSpikes(i) = size(temp,1);
    index_beg_spike = [index_beg_spike; temp];
    stemp = sprintf('N: %d Spike Detection Parameters\tSTD:\t%1.3f, M:\t%1.3f,\tTH:%d',nSpikes(i), fstd(i),fmean(i),std_threshold)
end

minTime = 15e-3; % 2ms
overlapWin = minTime/output.dt;
%% remove events that occur within minTime of eachother
index_beg_spike = sort(index_beg_spike);
index_beg_spike = index_beg_spike([1 (find(diff(index_beg_spike)>overlapWin)'+1)]);

clear temp;
if isempty(find(neighbors == -1))
    for i=1:size(neighbors,2)  %Extract neighbors (column number - offset in data) for triggerEE using ADC numbers
%         temp(i)  = find(output.nADCSamplingSeq == neighbors(triggerEE-nIntracellsignals,i));
        temp(i)  = find(output.nADCSamplingSeq == neighbors(1,i));
    end
    index_beg_spike = unique(index_beg_spike); %% Because there may be spikes detected by multiple electrodes
else
    temp = -1;    %Extract spikes
    %  it doesn't matter which triggerEE is used since spikes are from all
    %  electrodes
    index_beg_spike = unique(index_beg_spike); %% Because there may be spikes detected by multiple electrodes
    % note there is a potential problem because spikes within > 1 sample
    % are allowed
end

[allElectrodesXspikes spikeData temp skipped]  = xSpikesEE2(index_beg_spike,output.data,triggerEE,Ts,N_chn,maxpeakamp+fmean(find(EEgroup==triggerEE)),maxpeakwidth,temp,nIntracellsignals,WOI,fmean, fstd);
% [allElectrodesXspikes spikeData temp skipped]  = xSpikesEEgroup(EEgroup,index_beg_spike,output,nIntracellsignals,WOI);
for jj = 1:4
    Spikeskipped(jj) = sum(skipped == -1*jj);
end
clear skipped;
%% save data  to disk
if bsave

    dirtemp = '';
    figdesc = 'xSpikeWaveform';
    spath = [writedirheader dirtemp '\'];
    if (~isdir(spath))
        mkdir(spath);
    end
    %         temp = sprintf('%s%sEE%d',spath,figdesc,triggerEE);
    temp = sprintf('%s%sEE_%d',spath,figdesc,triggerEE);
    XspikeParams = struct('std',fstd,'mean',fmean,'std_threshold',std_threshold,'maxpeakamp',maxpeakamp,'maxpeakwidth',maxpeakwidth,'WOI',WOI);
    save(temp,'XspikeParams','neighbors','allElectrodesXspikes','spikeData','Spikeskipped','nSpikes','EEgroup');
end
% output.data = []; % savespace
% clear output.data;
% b_intra = 1;
%     %% find Intracellular Spikes
%     temp = int32(FindSpikes(output.data(:,indexoffset+triggerEE:N_chn:end), std_threshold*fstd(triggerEE) ,Ts,b_intra));
if bsortAmp
    %% parm = [P1 T1 P2 P1W T1W P2W T1F T1R P2R P2F];
    ele_param = xSpikeParamsCon(allElectrodesXspikes(:,:),output.dt,size(neighbors,2)+1)
  %% cell wherre each 
%     figure(15)
% %     plot(ele_param{1}(:,2),ele_param{2}(:,2),'.k');
%     plot(ele_param{1}(:,2)./ele_param{1}(:,3),ele_param{2}(:,2),'.k');
%     SCORE = [ele_param{1}(:,2),ele_param{2}(:,2)];
%AMP only
    parameters = single([ele_param{1}(:,2),ele_param{1}(:,3),ele_param{2}(:,2),ele_param{2}(:,3),ele_param{3}(:,2),ele_param{3}(:,3)]);    
 % AMP and width of spike1
    parameters = single([ele_param{1}(:,2),ele_param{1}(:,3),ele_param{2}(:,2),ele_param{2}(:,3),ele_param{3}(:,2),ele_param{3}(:,3),ele_param{1}(:,5),ele_param{1}(:,6)]);    
end
    % Overcluster data
    num_c = 26;
    n_units = 5;
    nPC = 20;

   s_sortbyPCA;
   %% ADD plot PCA 1:3
     SCORE = SCORE(1:end,1:nPC);

blinkage = 0;
if size(SCORE,2) > 2  %% can't do cluster linkage analysis if there are less then 3 dimensions
    blinkage = 1;
end
    % if kmeans_flag == 1
  s_kmeanskluster
breport = 1; 
%% Plot spikes from klusters
s_plotklusterspikes
% j = 10
% s_PSPtimes

% if blinkage
% %     s_plotclusterspikes
% else
%     s_plotmanualclusterspikes
% end
% 
% s_acceptUnits  % manually enter accepted klusters
% bfirst = 1;
% breport = 1;
%             bIPSC =0;
% s_findFirstSpikeInBurst %% extract 1st spike in burst
% % % % 
% xcorrAcceptedNeurons

%% look at the 10% that deviates the most from the mean

% 
% 


