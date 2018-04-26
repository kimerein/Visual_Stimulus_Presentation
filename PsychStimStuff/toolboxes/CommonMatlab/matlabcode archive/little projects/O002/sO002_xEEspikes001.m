% Extact spikes
% bin sweeps onto some size
% count spikes
close all;
% clear all;
total_extractedspikes = 0;
processedFiles = '';
Nspikes_in_File = [];
warning off 

% writedirheader = 'e:\My Documents\Scanziani Lab\Patch Data\';
writedirheader = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\mdata\';

%  writedirheader = readdirheader;
%  readdirheader =  'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\';
 readdirheader =  '\\Lindseysrig\My Documents\Patch Data\';
%  readdirheader = 'z:\Patch Data\';

%   readdirheader =  'D:\';
% readdirheader =  'e:\My Documents\Scanziani Lab\Patch Data\';
% readfileheader = ['ConRed2S3_2005_10_21_']
 readfileheader = ['Concatenatea_'];
% readfileheader = ['2005_10_21_'];
% savefileheader = ['2005_12_26_100-350-spikes'];
savefileheader = ['Concatenatea_']
readfilenumber  = [7];

for jj = [size(readfilenumber,2):-1:1]
    total_extractedspikes = 0;
    processedFiles = '';
    Nspikes_in_File = [];

    [filename] = createaxonfilename(readfileheader,readfilenumber(jj));
    pathfilename = strcat(readdirheader,filename);
    [sfilename] = creatematlabfilename(savefileheader,readfilenumber(jj));
    savefilename = strcat(writedirheader,sfilename);

    Chan = 1; %ONLY use 1 channel
%     inTs = 1/0.02e-3 ;%Sample rate
    inTs = 1/66e-6; %

    subsample = 1; %
    Ts = inTs/subsample; %subsampled
    sweeps = -1 ;%[450:700];
    %EXTRACT DATA
    tic

    %CHECK if same file is already loaded (saves time)
    b_loadnewdata =1;
    if(exist('lastfileloaded'))
        if (strcmp(pathfilename,lastfileloaded))
            b_loadnewdata = 0;
        end
    end
% b_loadnewdata =1
    if b_loadnewdata
%         sweeps = [300:400]
        % [my_data N_sweeps N_chan N_samples] = import_abf(pathfilename,sweeps,1/inTs)%;,Chan,subsample);
        pack;
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
        data =full(my_data(1:subsample:end,1+Chan+N_chn:N_chn:end)) ; %Data w/o
        warning('Skipped first sweep');% sperious values in beginning of first sweep.. so it is skipped
        N_sweepsExtracted = N_sweepsExtracted -1;
        %     warning('Inverted my_data');

        % SUBSET OF SWEEPS
        clear my_data; %save memory

        lastfileloaded = pathfilename;

    end %Loadnew data
    % pack;

    WOI = [1000e-6 1000e-6; 2000e-6 2000e-6];
    Waveform_size = [(WOI(1)+WOI(2))*Ts+1; (WOI(3)+WOI(4))*Ts+1];
    extracted_spikes = [];

    plot(data(:,1));
    hold all
    % plot(data(:,end),'-r');
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
pack;
index_beg_spike = FindSpikes(data, threshold,Ts,0);

    %% END of CHOOSE SPIKES
    toc;
    lastsweepnum = 0;
    tic%% EXTRACTDATA (%TODO funciton)
    if true % Runs faster w/o extracting spikes
        [tempXspikes spikeData sumspiketimeSw] = xSpikesEE(index_beg_spike,data,Ts,N_chn,maxpeakamp,maxpeakwidth);
        spiketime = spikeData(:,1);
        spiketimeSw = spikeData(:,[2:end]); % [sweepnum timeInsweep ISI]
        extracted_spikes = [extracted_spikes tempXspikes];
    else
        spiketime = index_beg_spike;
        extracted_spikes = [-1];
    end
    %filename must have 19 columns
    filename = [blanks(19-length(filename)) filename];
    processedFiles = [processedFiles; filename];
    %     Nspikes_in_File(ii) = size(index_beg_spike,1)
    Nextactedspikes = size(spiketime,1);
    toc

    % SAVE
    save(savefilename,...
        'processedFiles','extracted_spikes','Nextactedspikes','spiketime','N_sweeps','N_sweepsExtracted','N_samples','Ts','subsample','-mat');


    %%///////////////////////////////////////////////////

    %      Interval = [0 0];

%      run sO002_PLOT_print
%     pause;
end
