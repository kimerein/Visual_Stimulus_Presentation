% Extact spikes
% bin sweeps onto some size
% count spikes
close all;
% clear all;
total_extractedspikes = 0;
processedFiles = '';
Nspikes_in_File = [];

% writedirheader = 'e:\My Documents\Scanziani Lab\Patch Data\';
writedirheader = 'z:\Patch Data\';
%  writedirheader = readdirheader;
%  readdirheader =  'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\';
 readdirheader =  '\\Lindseysrig\My Documents\Patch Data\';
%  readdirheader = 'z:\Patch Data\';

%   readdirheader =  'D:\';
% readdirheader =  'e:\My Documents\Scanziani Lab\Patch Data\';
readfileheader = ['2006_01_19_'];
savefileheader = ['2006_01_19_'];
readfilenumber  = [4];
[filename] = createaxonfilename(readfileheader,readfilenumber(1));
pathfilename = strcat(readdirheader,filename);
[sfilename] = creatematlabfilename(savefileheader,readfilenumber(1));
savefilename = strcat(writedirheader,sfilename);

Chan = 4; %ONLY use 1 channel
% inTs = 1/0.02e-3 ;%Sample rate
inTs = 1/20e-6;
subsample = 1; %
Ts = inTs/subsample; %subsampled
sweeps = -1;
%EXTRACT DATA
tic

%CHECK if same file is already loaded (saves time)
b_loadnewdata =1;
if(exist('lastfileloaded'))
    if (strcmp(pathfilename,lastfileloaded))
        b_loadnewdata = 0;
    end
end

if b_loadnewdata
    % [my_data N_sweeps N_chan N_samples] = import_abf(pathfilename,sweeps,1/inTs)%;,Chan,subsample);
    [my_data N_sweeps N_chn N_samples] = import_abf(pathfilename,sweeps,1/inTs);%;,Chan,subsample);
    toc
    if sweeps ~=-1
        N_sweepsExtracted = size(sweeps,2);
    else
        N_sweepsExtracted = N_sweeps;
    end

    N_samples = N_samples/subsample;

    % N_chn = size(Chan,2);

    tic
    %Initialize Variables
    time_data = my_data(1:subsample:end,1);
    data =-my_data(1:subsample:end,1+Chan+N_chn:N_chn:end) ; %Data w/o
    warning('Skipped first sweep');% sperious values in beginning of first sweep.. so it is skipped
    N_sweepsExtracted = N_sweepsExtracted -1;
    warning('Inverted my_data');


    
    %IC data %first and last sweep
    temp = 1; %IC primary channel in files
    dataIC = zeros(size(my_data,1)/subsample,2);
    dataIC(:,1) = my_data(1:subsample:end,1+temp);
%     dataIC(:,2) = my_data(1:subsample:end,(N_sweepsExtracted-1)*N_chn + 1 +temp);
    %Current injection dat
    temp = 2; %IC secondary channel in files
    dataIC2 = zeros(size(my_data,1)/subsample,2);
    dataIC2(:,1) = my_data(1:subsample:end,1+temp);
    dataIC2(:,2) = my_data(1:subsample:end,(N_sweepsExtracted-1)*N_chn + 1 +temp);

    clear my_data; %save memory
    
    lastfileloaded = pathfilename;

end %Loadnew data
% pack;

WOI = [1000e-6 1000e-6; 2000e-6 2000e-6];
Waveform_size = [(WOI(1)+WOI(2))*Ts+1; (WOI(3)+WOI(4))*Ts+1];
extracted_spikes = [];

plot(data(:,1));
hold all
plot(data(:,end),'-r');
pause;
% USER ENTER Threshold
strArray = java_array('java.lang.String',2);
strArray(1) = java.lang.String('Enter threshold');
strArray(2) = java.lang.String('Enter maxpeakamp');
prompt = cell(strArray);
dlg_title = 'Event Detection Parameters'; num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines);%,defAns)
threshold = str2num(answer{1,1}); %10
maxpeakamp = str2num(answer{2,1});
maxpeakwidth =  1e-3;  % 
% FIND spikes
index_beg_spike = FindSpikes(data, threshold,Ts,0);

%% END of CHOOSE SPIKES
toc;
lastsweepnum = 0;
tic%% EXTRACTDATA (%TODO funciton)
[tempXspikes spikeData sumspiketimeSw] = xSpikesEE(index_beg_spike,data,Ts,N_chn,maxpeakamp,maxpeakwidth);
spiketime = spikeData(:,1);
spiketimeSw = spikeData(:,[2:end]); % [sweepnum timeInsweep ISI]

%filename must have 19 columns
filename = [blanks(19-length(filename)) filename];
processedFiles = [processedFiles; filename];
%     Nspikes_in_File(ii) = size(index_beg_spike,1)
extracted_spikes = [extracted_spikes tempXspikes];
Nextactedspikes = size(spiketime,1);
toc

% SAVE
save(savefilename,...
    'processedFiles','extracted_spikes','Nextactedspikes','spiketime','N_sweeps','N_sweepsExtracted','N_samples','Ts','subsample','-mat');


%%///////////////////////////////////////////////////

%Find INTERVAL time of Depolarization  %Secondary IC channel
% temp = dataIC2(500:end,1)+ abs(mean(dataIC2(500:550,1))) > .5;  % 500pA; %skip beginning because of garbage (unknown why) % offset to zero
toffset = round(size(dataIC2,1)/4);  %injected current only can occur in 2nd half of trace
temp = dataIC2(toffset:end,1);
temp1 = diff(temp);
tmax = find(temp1==max(temp1)); tmin = find(temp1==min(temp1));
if isempty(tmax) || isempty(tmin)
    Interval = [0 0];
% elseif (tmin(1)-tmax(1)) ~= sum(temp(tmax:tmin+1)>(temp(tmax+1)/2)) %current is continuously high between tmax and tmin
    Interval = [0 0];
else
    Interval = [tmax tmin]+toffset;
end
     Interval = [0 0];

run sO003_originalPlot
run sO003_PLOT_print
