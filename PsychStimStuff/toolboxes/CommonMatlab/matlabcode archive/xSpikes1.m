function [processedFiles,Nspikes_in_File,extracted_spikes] = xSpikes1(readfileheader,readfilenumbers)
%***********************************************************************
% Function
% [processedFiles,Nspikes_in_File,extracted_spikes] = xSpikes(readfileheader,readfilenumbers)
%
% Based on SpikeExtraction008
% Bassam Atallah
% Last change: 10/1/2005
% Reads .abf file,
%
% extracted_spikes(Amplitude,spike#,intra/header/extra)
% z = 1 intra, z= 3 extra
% Z = 2 intraheader, z = 4 extraheader
%
%***********************************************************************
% Initialize workspace
%Parameters for import_abf
xChan = [1 3]; %Channels to extract
Ts = 1/0.02e-3 ;%Sample rate 50kHz
%%%%
%%filename and path%%
readdirheader =  'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\';
% readdirheader =  'E:\My documents\Academic\Rotation\Scanziani\Patch Data\'; %HOME
%%%%

clc;
warning  off all;

%USE ONLY first time
processedFiles = char(zeros(1,19));
putative_type = [];
cell_number = [];
ID_SpikeLast = 0;
total_extractedspikes = 0;
extracted_spikes = [];

Nfiles_per_header = size(readfilenumbers,1); %modify for > 1 readfileheader

clear my_data;
i = 0;
Extra_Chn = [];
Intra_Chn = [];
nExtra_Chn = 0;
filename = '';
%//////////////////////////
for (j = 1: size(readfileheader,1))  %LOOP through headers
    filenumber_ind = i+1;
    for (i =filenumber_ind:filenumber_ind + Nfiles_per_header(j,1) -1)  % count through readfilenumbers for this readfileheader

        [filename] = createaxonfilename(readfileheader(j,:),readfilenumbers(i));
        filename = strcat(readdirheader,filename) ;

        Intra_Chn = 1;
        Extra_Chn = 2; %index in xChan
        nExtra_Chn = 1;
        N_chn = size(xChan,2);

        my_data_zind = 1; % DON"T KNOW why this variable exists
        clear my_data;
        [my_data(:,:,my_data_zind) N_sweeps] = import_abf(filename,-1,1/Ts,xChan);
        N_samples = size(my_data,1);

        %%% Plot abf files
        figure(figure_ID+i);
        for j = 1:size(xChan,2)
            subplot(2,ceil(size(xChan,2)/2),j);
            plot(my_Data(:,j+1:size(xChan,2):end),'r-');
            hold on;
            title(filename(end-18:end),'Interpreter','none');
        end


        % ---------------------------------------------------------------------
        % ----

        %PEAK extaction PARAMETERS
        %Window of interest around Beg_spike -500us + 2500us
        WOI = [1000e-6 1000e-6; 2000e-6 2000e-6];
        Waveform_size = [(WOI(1)+WOI(2))*Ts+1; (WOI(3)+WOI(4))*Ts+1];
        threshold = 10;
        minpeakwidth = .4e-3; %threshold to threshold
        maxpeakwidth = 2e-3;

        %Initialize Variables
        clear intraCell_data; clear extraCell_data; clear time_data
        time_data = my_data(:,1);
        intraCell_data =my_data(:,Intra_Chn+1:N_chn:end,my_data_zind) ;
        if (nExtra_Chn ~= 0)
            extraCell_data = my_data(:,Extra_Chn+1:N_chn:end,my_data_zind) ;
        else
            extraCell_data = zeros(size(my_data,1),Extra_Chn:N_chn:end,my_data_zind) ;
        end
        % Find spike from INTRAdata
        data_GT = intraCell_data > 10;
        temp = [diff(data_GT); zeros(1,size(data_GT,2))]; %insert a row because diff is small then the matrix it comes from
        index_beg_spike = find(temp == 1);
        index_end_spike = find(temp == -1);
        %% CHOOSE SPIKES with min and maxpeakwidth
        if (index_end_spike(1) < index_beg_spike(1))  % CASE Falling edge occurs before rising edge
            index_end_spike = index_end_spike(2:end);
        end
        if (size(index_end_spike(1),2) < size(index_beg_spike(1),2) ) %CASE Rising edge with no falling edge
            index_beg_spike = index_beg_spike(1:end-1);
        end
        peakwidth = index_end_spike - index_beg_spike;
        mask = (peakwidth > minpeakwidth*Ts) .* (peakwidth < maxpeakwidth*Ts); %Choose peaks with correct width %Indices that meet criteria have value 1
        index_beg_spike = nonzeros(mask .* index_beg_spike);
        %% END of CHOOSE SPIKES

        if ~isempty(index_beg_spike)

            % find Spike interval (beg to beg of Intracellular spike)
            index_beg_spike_sh = circshift(index_beg_spike,-1);
            index_beg_spike_sh(size(index_beg_spike_sh,1),1) = 0; %remove element shifted to the end of the vector
            interspike_interval = index_beg_spike_sh - index_beg_spike; %interval in index
            interspike_interval = circshift(interspike_interval,1);
            interspike_interval(1,1) = -1;

            %Note aligned temporally by 1st zero crossing of Intracellular spike
            for (ii  = 1:size(index_beg_spike,1)) % for all spikes found
                IndexES = total_extractedspikes + ii;
                Intra_cell = 1;
                Chn = Intra_Chn;
                z = 1;
                temp1 = intraCell_data (index_beg_spike(ii,1):index_beg_spike(ii,1)+ 1.6e-3*Ts)';
                %             Peak1 = find(intraCell_data (index_beg_spike(ii,1):index_beg_spike(ii,1)+ 80) == max(temp1)) + index_beg_spike(ii,1);
                Peak1 = find(temp1 == max(temp1)) + index_beg_spike(ii,1);
                %take first Peak
                tempXspikes(1:int16(Waveform_size(1)),IndexES,z) = intraCell_data(Peak1(1,1)-WOI(1)* Ts:Peak1(1,1)+WOI(2)* Ts)';
                Sweep =ceil(Peak1(1,1)/N_samples);
                if (Sweep > N_sweeps) %HACK not sure why it every can happen
                    Sweep = N_sweeps;
                end
                Timeofspike = Peak1(1,1) -N_samples*(Sweep -1);

                %EXTRACELLULAR SPIKE
                z = 3;
                tempXspikes(1:int16(Waveform_size(1)),IndexES,z) = extraCell_data (Peak1(1,1)-WOI(1)* Ts:Peak1(1,1)+WOI(2)* Ts)';

            end % for each spike
            %filename must have 19 columns
            filename = [blanks(19-length(filename)) filename];
            processedFiles = [processedFiles; filename];
            Nspikes_in_File(i) = size(index_beg_spike,1)
        end%ISEMPTY BEG_SPIKE
        extracted_spikes = [extracted_spikes tempXspikes];
    end
end
processedFiles
total_extractedspikes

%%Adjust baseline
% base = -60 - mean(extracted_spikes(1:10,:,1));
% Baseline = ones(size(extracted_spikes,1),1)*base;
% extracted_spikes(:,:,1) = extracted_spikes(:,:,1) + Baseline ;
%
% extracted_spikes(62,:,2) = base(1,:);  %Record factor in header

%%Scale Amplitude to 1
%{
temp = 1./max(extracted_spikes(:,:,1));
Scalar =  ones(size(extracted_spikes,1),1)*temp;
extracted_spikes(:,:,1) = extracted_spikes(:,:,1).*Scalar ;
%}
% scal = -1./min(extracted_spikes(:,:,3));
% Scalar =  ones(size(extracted_spikes,1),1)*scal;
% extracted_spikes(:,:,3) = extracted_spikes(:,:,3).*Scalar ;
% extracted_spikes(63,:,4) = scal(1,:); %Record factor in header

