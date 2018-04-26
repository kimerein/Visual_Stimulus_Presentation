% clear all
defineDir
warning off
expttype = 'KinateOsc';
writedirpath = [DATAANAL_DIR expttype '\'];
writedirheader = [DATAANAL_DIR expttype '\'];
breport = 1;
bsave = 1;
savetype = 'emf'
exptnum = 's1';
% exptnum = 's2c1';
% exptnum = '2s6c1';
% xcolor = '-g'; v = 7.5; % IPSCs
% xcolor = '-r'; v = -82; % EPSCs
xcolor = '-c'; v = -62;% mixed
bVC = 0;
% output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_07_20_0011.abf',-1,1);
output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_08_08_0006.abf',-1,1);

if ~strcmp(lastexptnum,exptnum)
close all
end



thres = -0.01; %% MANUALLY DEFINE threshold for oscilation detection
 LPFilter = 60 ;%% low pass filter of LFP
 ISIrange = [1/20 1/50]*1000; % ms
odata = struct();
odata.expttype = expttype;
odata.exptnum = [exptnum '_' output.sfilename(1:strfind(output.sfilename,'.')-1)]; %% UNIQUE KEY
odata.sfilename = output.sfilename;
odata.dt = output.dt;
odata.holdVoltage = v;
odata.VC = bVC;
odata.thres = thres;
odata.LPFilter = LPFilter;
odata.ISIrange = ISIrange;
odata.analysisDate = datestr(now);

fileindex =  output.sfilename(max(strfind(output.sfilename,'_'))+1:max(strfind(output.sfilename,'_'))+4);
writedirheader = [writedirheader output.sfilename(1:max(strfind(output.sfilename,'_'))-1) '_' exptnum '\'];
if ~isdir(writedirheader)
    mkdir(writedirheader);
end;
indexoffset = 1; b_intra= 0;

  
%%LFP
lfpdata = double(output.data(:,indexoffset+2:output.Nchan:end))*output.gain(2);
% lfpdata = lfpdata/mean(std(lfpdata));
lfpdata = (lfpdata-mean(mean(lfpdata)));
intradata = double(output.data(:,indexoffset+1:output.Nchan:end));
% intradata = intradata - mean(mean(intradata));
% intradata = intradata/min(std(intradata));

%% low pass filter LFP data
      [B,A] = butter(2,2*LPFilter*output.dt,'low');
prodata =   filtfilt(B,A,lfpdata);
bdebug = 0;
if bdebug
 plot(lfpdata(1:5000,2))
 hold on
 plot(t(1:5000,2),'r')
end
%% find trough of each oscilation
  %% reject events out side of ISIrange
  %% extract WOI in Intradata after event
findOSC_detectOSC

osc = bin2(Strough',1/output.dt/1000);

%%% spike times  ( from s_tempSpikeSort2)
nIntracellsignals =0; % #ADC channels which were Intracellular (Code assumes they are contiguous from 1 to N)
indexoffset = 1  ; % #Offset to ADC channels cause data starts with time column
neighbors = [1] % ADC number electrodes that are used with other electrode in spikesorting (waveforms will be concatinated)
EEgroup = [1];


%% CREATE DIR for data
writedirheader = [writedirheader output.sfilename(1:end-4) '\' 'EE' strrep(num2str(EEgroup),'  ','_') '\'];
if ~isdir(writedirheader)
    mkdir(writedirheader);
end;

figheader = '';
%%% EXTRACTED from DATA
Ts = 1/output.dt; % Sampling Frequency (Typically 20kHz)
N_chn = output.Nchan; % Number of ADC Channels recorded from
N_samples = output.Nsamples; % Number of ADC Channels recorded from
Nsweeps = output.Nsweeps;
sourcefile = output.sfilename;
sw = max(output.xSweeps); % sweeps to set max time
sourcefile = output.sfilename;
for i=1:nIntracellsignals
    IntraCellData = output.data(:,1+i:output.Nchan:end);
end

%%% Spike Extraction Parameters
WOI = [800e-6 800e-6; 1100e-6 1100e-6];%Window of interest around Beg_spike -1000us + 2000us
std_threshold =3; %number of std for threshold
b_intra = 0; % trigger on negative going spike
maxpeakamp = 10000*std_threshold;  %in std
maxpeakwidth =  2e-3;


index_beg_spike = 0; temp1 = []; temp2 = []; temp2 = []; temp = []; fstd = []; fmean = [];
% for triggerEE=nIntracellsignals+1:N_chn  % all Extracellular signals , threshold
for i = 1: size(EEgroup,2)  % all Extracellular signals in 1 tetrode
    triggerEE = EEgroup(i);
    if ~isinteger(output.data)
        fstd(i) =  mean(std(output.data(:,indexoffset+triggerEE+N_chn:N_chn:5*N_chn)));
    else
        fstd(i) =  mean(int16std(output.data(:,indexoffset+triggerEE+N_chn:N_chn:5*N_chn)));
    end

    fmean(i) =  mean(mean(output.data(:,indexoffset+triggerEE+N_chn:N_chn:5*N_chn)));

    %%% CHANGE code so spikes are in std units
    tic
    %% REmember threshold is mean - std*std_threshold
    temp = [temp; int32(FindSpikes2(output.data(:,indexoffset+triggerEE:N_chn:end), (fmean(i)-(std_threshold*fstd(i))) ,Ts,b_intra));];
    toc
    nSpikes(i) = size(temp,1);
    index_beg_spike = [index_beg_spike; temp];
    stemp = sprintf('N: %d Spike Detection Parameters\tSTD:\t%1.3f, M:\t%1.3f,\tTH:%d',nSpikes(i), fstd(i),fmean(i),std_threshold)
end

minTime = 1e-3; % 2ms
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


%%%%%%%%%%%%%%%555555
spk = bin2(index_beg_spike,1/output.dt/1000);% bined to 1ms
a = xcorr(double(spk), double(osc),int32(200)); % bined to 1ms
plot(([1:size(a,2)]-ceil(size(a,2)/2)),a)

