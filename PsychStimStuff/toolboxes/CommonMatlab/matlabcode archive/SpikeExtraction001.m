%***********************************************************************
% SpikeExtraction.m
%
% Bassam Atallah
% Last change: 8/01/2005
% Reads .atf file, Extacts header and spikes.
% Spikes placed in
%                 extracted_spikes(Amplitude,spike#,intra/header/extra/header)
% z = 1 intra, z= 3 extra
% Z = 2 intraheader, z = 4 extraheader
%
%***********************************************************************




% Initialize workspace

clc;
warning  off all;

%Initialize  Plots

%close 4

%REad LastSpike
%KEY database parameters

%USE ONLY first time
ID_SpikeLast = 0;
total_extractedspikes = 0;
clear extracted_spikes;
%load('E:\My documents\Academic\Rotation\Scanziani\Data Analysis\072005 Clustering\extracted_spikes.mat');
%processedFiles
%total_extractedspikes
%describe datastucture of extacted_spikes

%INPUT File parameters
dirheader =  'e:\My Documents\Academic\Rotation\Scanziani\Ascii\';
fileheader = ['2005_07_12_';'2005_07_14_']
Nfiles_per_header = [3 ; 2]
PROTOCOL = 1; %ICP = 1, ICR = 2
filenumbers = [4 11 18 47 56]; %file index for day
CELL_NUM = [07120501 07120502 07120503 07140505 07140506] %dateCell# 
PUT_TYPE = [1 2 1 1 0]% 1= PYR, 2 = Basket 0 = Other

MAX_CHN = 16; %max ADC channels in datafile header
HEADER_LENGTH = 50;
%NUMBER_SPIKES = 250;
if(sum(Nfiles_per_header) ~= size(filenumbers,2))
    Warning = 'CHECK INPUT File parameters   Nfiles_per_header ~= size(filenumbers,2)'
end
clear my_data;
i = 0;
%//////////////////////////
for (j = 1: size(fileheader,1))
   filenumber_ind = i+1;
    for (i =filenumber_ind:filenumber_ind + Nfiles_per_header(j,1) -1)  % count through filenumbers for this fileheader
        %Set Filename
        if (filenumbers(i)<1000)
            filename = strcat(fileheader(j,:),'0');
        end
        if (filenumbers(i)<100)
            filename = strcat(filename,'0');
        end
        if (filenumbers(i)<10)
            filename = strcat(filename,'0');
        end
        filename = sprintf('%s%d.atf',filename,filenumbers(i))
        file = strcat(dirheader,filename) ;
        fid = fopen(file,'r');

        % Load header
        %READING HEADER
        temp = textscan(fid,'%*f %f\n','headerLines',1);
        N_col = cell2mat(temp);
        temp = textscan(fid,'"SignalsExported=','headerLines',7);
        temp = textscan(fid,'%s""'); %if I don't make it into a string the cell have '' rather then are empty... (dont' know why)
        temp = textscan(cell2mat(temp{1,1}),'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
        N_chn = MAX_CHN - sum(cellfun('isempty',temp));
        Chn_name = char(zeros(N_chn,6));
        for ii = 1: N_chn
            a = cell2mat(temp{1,ii});
            Chn_name(ii,1:6) = a(1,1:6)  %how do I copy strings
            if ( isempty(regexp(Chn_name(ii,:),'IC')) ~= 1)
                Intra_Chn = ii
            end
            if ( isempty(regexp(Chn_name(ii,:),'ff')) ~= 1)
                Extra_Chn = ii
            end
        end
        N_sweeps = (N_col-1)/N_chn
        %need to go back to the beginning of the file
        fid = fopen(file,'r');
        temp = textscan(fid,'"SweepStartTimesMS=%s\n','headerLines',7);
        str_t = cell2mat(temp{1,1});
        temp = textscan(str_t,'%f,');
        Sweep_Times = cell2mat(temp); %row mat
        % END READ HEADER

        % LOAD data
        %toggle saving all data
        my_data_zind = 1; % i;
        clear my_data;
        
        my_data(:,:,my_data_zind) = dlmread(file,'\t',12,0);
        N_samples = size(my_data,1);
        Ts = 1/( my_data(2,1,my_data_zind) - my_data(1,1,my_data_zind)); %usually 50kHz


        % ---------------------------------------------------------------------
        % ----

        %PEAK extaction PARAMETERS
        %Window of interest around Beg_spike -500us + 2500us
        WOI = [1000e-6 1000e-6; 2000e-6 2000e-6];
        Waveform_size = [(WOI(1)+WOI(2))*Ts+1; (WOI(3)+WOI(4))*Ts+1];
        threshold = 10;
        %minpeakwidth = 2000e-6;
        %maxpeakwidth = 4000e-6

        %Initialize Variables
        clear intraCell_data; clear extraCell_data; clear time_data
        time_data = my_data(:,1);
        intraCell_data =my_data(:,Intra_Chn+1:N_chn:end,my_data_zind) ;
        extraCell_data = my_data(:,Extra_Chn+1:N_chn:end,my_data_zind) ;
    
        %fix incorrect GAIN setting
        extraCell_data = extraCell_data * 10;

        % Find spike from INTRAdata
        data_GT = intraCell_data > 10;
        temp = [diff(data_GT); zeros(1,size(data_GT,2))]; %insert a row because diff is small then the matrix it comes from
        index_beg_spike = find(temp == 1);

        % find Spike interval (beg to beg of Intracellular spike)
        index_beg_spike_sh = circshift(index_beg_spike,-1);
        index_beg_spike_sh(size(index_beg_spike_sh,1),1) = 0; %remove element shifted to the end of the vector
        interspike_interval = index_beg_spike_sh - index_beg_spike; %interval in index
        interspike_interval = circshift(interspike_interval,1);
        interspike_interval(1,1) = -1;

        %Note aligned temporally by 1st zero crossing of Intracellular spike
        for (ii  = 1:size(index_beg_spike,1)) % for all spikes found
            IndexES = total_extractedspikes + ii;
            for(iii = 1:2) %deal with intra/extra independently
                if(iii == 1)
                    Intra_cell = 1;
                    Chn = Intra_Chn;
                    z = 1;
                    temp1 = intraCell_data (index_beg_spike(ii,1):index_beg_spike(ii,1)+ 80)';
                    Peak1 = find(intraCell_data (index_beg_spike(ii,1):index_beg_spike(ii,1)+ 80) == max(temp1)) + index_beg_spike(ii,1);
                    %take first Peak
                    extracted_spikes(1:int16(Waveform_size(1)),IndexES,z) = intraCell_data(Peak1(1,1)-WOI(1)* Ts:Peak1(1,1)+WOI(2)* Ts)';
                    Sweep =ceil(Peak1(1,1)/N_samples);
                    Timeofspike = Peak1(1,1) -N_samples*(Sweep -1);
                else
                    Intra_cell = 0;
                    Chn = Extra_Chn;
                    z = 3;
                    temp2 = extraCell_data (index_beg_spike(ii,1)-50:index_beg_spike(ii,1)+ 30)';
                    Trough1 = find(extraCell_data (index_beg_spike(ii,1)-50:index_beg_spike(ii,1)+ 30) == min(temp2)) + index_beg_spike(ii,1)-50;
                    %take first Trough
                    extracted_spikes(1: int16(Waveform_size(1)),IndexES,z) = extraCell_data (Trough1(1,1)-WOI(1)* Ts:Trough1(1,1)+WOI(2)* Ts)';
                    Timeofspike(1,1) = Trough1(1,1) -N_samples*(Sweep -1);
                end
                %DO NOT rearrange the order
                % VER 001
                ID_SpikeLast = ID_SpikeLast + 1;
                extracted_spikes(1,IndexES,z+1) = HEADER_LENGTH;
                extracted_spikes(2,IndexES,z+1) = ID_SpikeLast ;
                extracted_spikes(3,IndexES,z+1) = CELL_NUM(1,i);
                extracted_spikes(4,IndexES,z+1) = filenumbers(1,i);%* filename data_XXXX  
                extracted_spikes(5,IndexES,z+1) = ii;%*spike number in file
                extracted_spikes(6,IndexES,z+1) = size(index_beg_spike,1);%*total spikes in file
                
          %{
                extracted_spikes(5,IndexES,z+1) = PreAMP; 
                extracted_spikes(5,IndexES,z+1) = AMP;                 
                extracted_spikes(5,IndexES,z+1) = Stim (yes/no); 
                extracted_spikes(5,IndexES,z+1) = (Stim electrode)       
                extracted_spikes(5,IndexES,z+1) = (Stim type bi phasic voltage)
                extracted_spikes(5,IndexES,z+1) = (Stim amplitude)
                extracted_spikes(5,IndexES,z+1) = Recording Electrodes
                extracted_spikes(5,IndexES,z+1) = HiPASS; 
                extracted_spikes(5,IndexES,z+1) = LoPASS; 
                extracted_spikes(5,IndexES,z+1) = Gain;
                %}
                
                extracted_spikes(20,IndexES,z+1) = PROTOCOL;
                extracted_spikes(30,IndexES,z+1) = Ts;
                extracted_spikes(31,IndexES,z+1) = N_samples;
                extracted_spikes(32,IndexES,z+1) = Sweep;
                extracted_spikes(33,IndexES,z+1) = Sweep_Times(Sweep); %time in sec from beginning of experiment
           
                extracted_spikes(41,IndexES,z+1) = Chn;
                extracted_spikes(42,IndexES,z+1) = Intra_cell;
                extracted_spikes(43,IndexES,z+1) = PUT_TYPE(1,i);
             
                extracted_spikes(60,IndexES,z+1) = interspike_interval(ii,1);
                extracted_spikes(61,IndexES,z+1) = Timeofspike;

                %Baseline Adjust (SEE below
                %Scalar
            end

        end % for each spike
        total_extractedspikes = total_extractedspikes + ii;

     
        processedFiles(size(processedFiles,1)+1,:) = filename;
        Nspikes_in_File(i) = size(index_beg_spike,1)


    end
end
processedFiles
total_extractedspikes

%Adjust baseline
base = -60 - mean(extracted_spikes(1:10,:,1));
Baseline = ones(size(extracted_spikes,1),1)*base;
extracted_spikes(:,:,1) = extracted_spikes(:,:,1) + Baseline ;

extracted_spikes(62,:,2) = base(1,:);  %Record factor in header

%Scale Amplitude to 1
%{
temp = 1./max(extracted_spikes(:,:,1));
Scalar =  ones(size(extracted_spikes,1),1)*temp;
extracted_spikes(:,:,1) = extracted_spikes(:,:,1).*Scalar ;
%}
scal = -1./min(extracted_spikes(:,:,3));
Scalar =  ones(size(extracted_spikes,1),1)*scal;
extracted_spikes(:,:,3) = extracted_spikes(:,:,3).*Scalar ;
extracted_spikes(63,:,4) = scal(1,:); %Record factor in header

save('E:\My documents\Academic\Rotation\Scanziani\Data Analysis\072005 Clustering\extracted_spikes.mat',...
    'extracted_spikes', 'total_extractedspikes', 'ID_SpikeLast', 'processedFiles', 'Nspikes_in_File', '-mat');
%{
Example File:

ATF	1.0
8	101
"AcquisitionMode=Episodic Stimulation"
"Comment="
"YTop=100,1"
"YBottom=-100,-1"
"SyncTimeUnits=1000"
"SweepStartTimesMS=2625.000,5625.000,8623.000,11625.000,14625.000,17624.000,20624.000,23620.000,26627.000,29621.000,32626.000,35626.000,38757.000,41623.000,44626.000,47621.000,50627.000,53628.000,56627.000,59627.000,62626.000,65627.000,68623.000,71627.000,74630.000,77629.000,80628.000,83625.000,86625.000,89627.000,92627.000,95630.000,98630.000,101630.000,104624.000,107629.000,110625.000,113662.000,116660.000,119659.000,122663.000,125662.000,128655.000,131656.000,134660.000,137662.000,140661.000,143656.000,146661.000,149662.000"
"SignalsExported=BA_IClamp,ff_extra_"
"Signals="	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"	"BA_IClamp"	"ff_extra_"
"Time (s)"	"Trace #1 (mV)"	"Trace #1 (mV)"	"Trace #2 (mV)"	"Trace #2 (mV)"	"Trace #3 (mV)"	"Trace #3 (mV)"	"Trace #4 (mV)"	"Trace #4 (mV)"	"Trace #5 (mV)"	"Trace #5 (mV)"	"Trace #6 (mV)"	"Trace #6 (mV)"	"Trace #7 (mV)"	"Trace #7 (mV)"	"Trace #8 (mV)"	"Trace #8 (mV)"	"Trace #9 (mV)"	"Trace #9 (mV)"	"Trace #10 (mV)"	"Trace #10 (mV)"	"Trace #11 (mV)"	"Trace #11 (mV)"	"Trace #12 (mV)"	"Trace #12 (mV)"	"Trace #13 (mV)"	"Trace #13 (mV)"	"Trace #14 (mV)"	"Trace #14 (mV)"	"Trace #15 (mV)"	"Trace #15 (mV)"	"Trace #16 (mV)"	"Trace #16 (mV)"	"Trace #17 (mV)"	"Trace #17 (mV)"	"Trace #18 (mV)"	"Trace #18 (mV)"	"Trace #19 (mV)"	"Trace #19 (mV)"	"Trace #20 (mV)"	"Trace #20 (mV)"	"Trace #21 (mV)"	"Trace #21 (mV)"	"Trace #22 (mV)"	"Trace #22 (mV)"	"Trace #23 (mV)"	"Trace #23 (mV)"	"Trace #24 (mV)"	"Trace #24 (mV)"	"Trace #25 (mV)"	"Trace #25 (mV)"	"Trace #26 (mV)"	"Trace #26 (mV)"	"Trace #27 (mV)"	"Trace #27 (mV)"	"Trace #28 (mV)"	"Trace #28 (mV)"	"Trace #29 (mV)"	"Trace #29 (mV)"	"Trace #30 (mV)"	"Trace #30 (mV)"	"Trace #31 (mV)"	"Trace #31 (mV)"	"Trace #32 (mV)"	"Trace #32 (mV)"	"Trace #33 (mV)"	"Trace #33 (mV)"	"Trace #34 (mV)"	"Trace #34 (mV)"	"Trace #35 (mV)"	"Trace #35 (mV)"	"Trace #36 (mV)"	"Trace #36 (mV)"	"Trace #37 (mV)"	"Trace #37 (mV)"	"Trace #38 (mV)"	"Trace #38 (mV)"	"Trace #39 (mV)"	"Trace #39 (mV)"	"Trace #40 (mV)"	"Trace #40 (mV)"	"Trace #41 (mV)"	"Trace #41 (mV)"	"Trace #42 (mV)"	"Trace #42 (mV)"	"Trace #43 (mV)"	"Trace #43 (mV)"	"Trace #44 (mV)"	"Trace #44 (mV)"	"Trace #45 (mV)"	"Trace #45 (mV)"	"Trace #46 (mV)"	"Trace #46 (mV)"	"Trace #47 (mV)"	"Trace #47 (mV)"	"Trace #48 (mV)"	"Trace #48 (mV)"	"Trace #49 (mV)"	"Trace #49 (mV)"	"Trace #50 (mV)"	"Trace #50 (mV)"
0	-66.5253	-6.40869e-4	-66.5253	-3.05176e-4	-66.6077	-7.32422e-4	-66.5588	2.13623e-4	-66.748	-0.00149536	-66.7633	4.27246e-4	-67.1173	2.13623e-4	-67.0746	-0.00115967	-66.8121	-3.96729e-4	-66.7236	-8.23975e-4	-66.9556	-7.32422e-4	-66.983	-3.96729e-4	-67.4927	-8.23975e-4	-67.5903	-8.54492e-4	-67.6025	-4.88281e-4	-67.2516	-4.88281e-4	-67.2333	-1.2207e-4	-66.8426	8.54492e-4	-67.2089	-6.10352e-4	-67.511	1.2207e-4	-67.279	-0.00119019	-67.3828	-3.66211e-4	-67.0685	4.27246e-4	-67.099	-7.62939e-4	-66.6931	-2.44141e-4	-66.5924	1.52588e-4	-66.8457	2.13623e-4	-66.3666	-2.44141e-4	-66.6504	-7.32422e-4	-66.9159	-6.10352e-5	-66.7419	-6.10352e-5	-66.8457	-0.00106812	-66.7267	2.44141e-4	-66.2659	-6.10352e-4	-66.5619	-8.23975e-4	-66.861	-3.96729e-4	-67.0929	-3.66211e-4	-66.8793	-7.01904e-4	-66.9556	-6.71387e-4	-66.925	-2.13623e-4	-66.9006	-6.40869e-4	-66.8915	5.49316e-4	-67.0746	-3.05176e-4	-66.4398	-5.18799e-4	-66.4368	-0.00112915	-66.3208	-0.00125122	-66.1957	3.35693e-4	-66.0644	-9.15527e-4	-65.9943	-3.66211e-4	-65.8813	3.05176e-5
2e-5	-66.5314	-3.96729e-4	-66.5131	-6.40869e-4	-66.6199	-6.40869e-4	-66.5619	-2.74658e-4	-66.745	-0.00125122	-66.7602	-7.62939e-4	-67.1112	-0.00140381	-67.0593	-0.0010376	-66.8091	-0.00131226	-66.7297	-5.49316e-4	-66.9403	3.05176e-4	-66.9861	-0.00146484	-67.4805	1.83105e-4	-67.5873	-6.10352e-4	-67.5964	-7.62939e-4	-67.2577	-5.49316e-4	-67.2333	1.2207e-4	-66.8304	-9.46045e-4	-67.2089	-6.71387e-4	-67.5201	-0.00112915	-67.2668	-8.54492e-4	-67.3767	-8.23975e-4	-67.0563	1.52588e-4	-67.1081	-3.05176e-5	-66.6992	-5.49316e-4	-66.5893	-7.62939e-4	-66.8487	-0.00158691	-66.3635	-2.44141e-4	-66.6382	-9.15527e-4	-66.9098	-2.74658e-4	-66.7358	-6.71387e-4	-66.8487	3.05176e-4	-66.7297	3.05176e-4	-66.2659	-5.79834e-4	-66.5558	-1.2207e-4	-66.8671	-1.52588e-4	-67.0807	-6.10352e-5	-66.8732	-0.00112915	-66.9495	-8.8501e-4	-66.9098	5.49316e-4	-66.8915	-5.49316e-4	-66.8915	-6.40869e-4	-67.0654	3.66211e-4	-66.4429	-3.35693e-4	-66.4307	-0.00143433	-66.3208	-9.15527e-4	-66.2018	-3.96729e-4	-66.0614	6.10352e-5	-65.9943	-7.32422e-4	-65.8844	-5.79834e-4
%}

%{
Extract these parameters from Header:
II
# Traces
Channels (name and number)
sample rate
%}
%}