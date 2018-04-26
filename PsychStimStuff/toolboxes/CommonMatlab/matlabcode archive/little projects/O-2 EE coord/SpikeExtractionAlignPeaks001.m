%***********************************************************************
% Based on SpikeExtraction008
%
% Bassam Atallah
% Last change: 9/29/2005
% Reads .abf file, 
% xSpikes placed in
%                 extracted_spikes(Amplitude,spike#,intra/header/extra/header)
% z = 1 intra, z= 3 extra
% Z = 2 intraheader, z = 4 extraheader
%
%***********************************************************************
% Initialize workspace

clc;
warning  off all;

%REad LastSpike
%KEY database parameters

% savefilename = 'E:\My documents\Academic\Rotation\Scanziani\Data Analysis\ExtraIntra\extracted_spikes.mat'
savefileheader = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\mdata\';
savefilename = '270905_B01.mat';
savefilename = strcat(savefileheader,savefilename);
%USE ONLY first time
processedFiles = char(zeros(1,19));
putative_type = [];
cell_number = [];
ID_SpikeLast = 0;
total_extractedspikes = 0;
extracted_spikes = [];
%%%%%%%%%%%%%%%%%%%%%%%%
%load('E:\My documents\Academic\Rotation\Scanziani\Data Analysis\072005 Clustering\extracted_spikes.mat');
%processedFiles
% total_extractedspikes
%describe datastucture of extacted_spikes

%INPUT File parameters
dirheader =  'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\'; 
% dirheader =  'E:\My documents\Academic\Rotation\Scanziani\Patch Data\'; %HOME
fileheader = ['2005_09_28_'];
% filenumbers = [18:22 24:30]; %file index for day
% Nfiles_per_header = [5];
%for one file header
Nfiles_per_header = size(filenumbers,1);

CONFIG = 2; %Recording configuration (electrodes, amplifiers, filtering) SEE e:\My Documents\Academic\Rotation\Scanziani\Data Analysis\Cell Summary (version 1).xls
PROTOCOL = 3; %ICP = 1, ICR = 2, ICP2 = 3
CELL_NUM = [05921000]; %dateCell# 
PUT_TYPE = [3];% 1= PYR, 2 = Regular spiking 3 = Fast Spiking 0 = Other
filename = '';
MAX_CHN = 16; %max ADC channels in datafile header
HEADER_VERSION = 50;

if(sum(Nfiles_per_header) ~= size(filenumbers,2) & sum(Nfiles_per_header) ~= size(filenumbers,1) & sum(Nfiles_per_header) ~= size(CELL_NUM) & sum(Nfiles_per_header) ~= size(PUT_TYP))
    Warning = 'CHECK INPUT File parameters   Nfiles_per_header ~= size(filenumbers,2)'
end
clear my_data;
i = 0;
Extra_Chn = [];
Intra_Chn = [];
nExtra_Chn = 0;
%//////////////////////////
for (j = 1: size(fileheader,1))
   filenumber_ind = i+1;
    for (i =filenumber_ind:filenumber_ind + Nfiles_per_header(j,1) -1)  % count through filenumbers for this fileheader
        %Set Filename
        if ~isempty(regexp(fileheader(j,:),'_'))
           if (filenumbers(i)<1000)
            filename = strcat(fileheader(j,:),'0');
           end
        
        else
          filename = fileheader(j,:);

        end
        if (filenumbers(i)<100)
            filename = strcat(filename,'0');
        end
        if (filenumbers(i)<10)
            filename = strcat(filename,'0');
        end
        filename = sprintf('%s%d.abf',filename,filenumbers(i))
        file = strcat(dirheader,filename) ;
        Intra_Chn = 1;
        Extra_Chn = 2; %index in xChan
        nExtra_Chn = 1;
        xChan = [1 3]; %Channels to extract
        N_chn = size(xChan,2);
        Ts = 1/0.02e-3 ;%Sample rate 50kHz
        my_data_zind = 1; % DON"T KNOW why this variable exists
        clear my_data;
        [my_data(:,:,my_data_zind) N_sweeps] = import_abf(file,-1,1/Ts,xChan);
        N_samples = size(my_data,1);
        
        
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
        mask = (peakwidth > minpeakwidth*Ts) .* (peakwidth < maxpeakwidth*Ts); %Indices that meet criteria have value 1
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

%Adjust baseline
% base = -60 - mean(extracted_spikes(1:10,:,1));
% Baseline = ones(size(extracted_spikes,1),1)*base;
% extracted_spikes(:,:,1) = extracted_spikes(:,:,1) + Baseline ;
% 
% extracted_spikes(62,:,2) = base(1,:);  %Record factor in header

%Scale Amplitude to 1
%{
temp = 1./max(extracted_spikes(:,:,1));
Scalar =  ones(size(extracted_spikes,1),1)*temp;
extracted_spikes(:,:,1) = extracted_spikes(:,:,1).*Scalar ;
%}
% scal = -1./min(extracted_spikes(:,:,3));
% Scalar =  ones(size(extracted_spikes,1),1)*scal;
% extracted_spikes(:,:,3) = extracted_spikes(:,:,3).*Scalar ;
% extracted_spikes(63,:,4) = scal(1,:); %Record factor in header

save(savefilename,'processedFiles','Nspikes_in_File','extracted_spikes','EE_coord', '-mat');

%save mean with name

%****OBJECTIVE online analysis.
