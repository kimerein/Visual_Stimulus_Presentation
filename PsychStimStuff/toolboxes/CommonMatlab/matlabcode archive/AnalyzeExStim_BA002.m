%***********************************************************************
% AnalyzeExStim.m
%
% Flavio Frohlich
% Last change: 7/10/2005
% Reads .atf file, basic analysis of extracellular traces with activity
% evoked by extracellular stimulation.
%
%***********************************************************************




% Initialize workspace

close all;
clc;

my_spikes = [];
candidates_proc = [];
my_ind = [];
HEADER_LENGTH = 50;
NUMBER_SPIKES = 250;
dirheader =  'e:\My Documents\Academic\Rotation\Scanziani\Ascii\';
fileheader = '2005_07_12_';
filenumbers = [18];

for (i =1: size(filenumbers,2))
    %Set Filename
    file = strcat(dirheader,fileheader) ;
    if (filenumbers(i)<1000)
        file = strcat(file,'0');
    end
    if (filenumbers(i)<100)
        file = strcat(file,'0');
    end
    if (filenumbers(i)<10)
        file = strcat(file,'0');
    end
    file = sprintf('%s%d.atf',file,filenumbers(i))
    % Load header
    %textread(file, 
    % Load data
    %my_data = dlmread(file,'\t',12,0);
    fid = fopen(file);

    fscanf(fid,
    [N_samples N_sweeps] = size(my_data);

    Ts = 1/( my_data(2,1) - my_data(1,1)); %usually 50kHz
    %Window of interest around Beg_spike -500us + 2500us
    WOI = [1000e-6 1000e-6; 2500e-6 2500e-6];
    Waveform_size = [(WOI(1)+WOI(2))*Ts+1; (WOI(3)+WOI(4))*Ts+1];

    % A = my_data(index_beg_conf_spike_1(:,1)-WOI(1)* Ts:index_beg_conf_spike_1(:,1)+WOI(2)* Ts,col)
    extracted_spikes_Intra = zeros(int16(Waveform_size(1)),NUMBER_SPIKES + HEADER_LENGTH,size(filenumbers,2));
    extracted_spikes_Extra = zeros(int16(Waveform_size(2)),NUMBER_SPIKES + HEADER_LENGTH,size(filenumbers,2));
%Datastructure:
% 50 rows header
    
%Database Structure
%Waveform col
% HEADER_LENGTH
% Filename
% Protocol (Pulse Burst)
% Channel
% Sweep

% Time in Sweep
% Ts
% Interspike interval
% Baseline Adjustment
% Scale factor

%File Array
% Filename
% Protocol
% Number of Channels
% Channel Names
% # Sweep
% Times of Sweeps
% Ts
    % -------------------------------------------------------------------------

    %
    % Find Peaks (*Need to think of a way w/o a for loop)


    threshold = 10;
    minpeakwidth = 2000e-6;
    maxpeakwidth = 4000e-6
    indexIncrementminpeakwidtth = floor(minpeakwidth * Ts); % converted ms into indicies
    mii = 1;
    max_index = zeros(10);
    col = 2;
    %Initialize Variables
    clear intraCell_data;
    clear extraCell_data;
    time_data = my_data(:,1);
    intraCell_data = my_data(:,2:2:end);
    extraCell_data = my_data(:,3:2:end);
    %clear my_data;

    % Find spike
    data_GT = intraCell_data > 10;
    temp = [diff(data_GT); zeros(1,size(data_GT,2))]; %insert a row because diff is small then the matrix it comes from
    index_beg_spike = find(temp == 1);

    % find Spike interval (beg to beg of Intracellular spike)
    index_beg_spike_sh = circshift(index_beg_spike,-1);
    index_beg_spike_sh(size(index_beg_spike_sh,1),1) = 0; %remove element shifted to the end of the vector
    spike_interval = index_beg_spike_sh - index_beg_spike; %interval in index
    spike_interval(size(spike_interval,1)) = -1;


    %Note aligned temporally by 1st zero crossing of Intracellular spike
    for (ii  =1:size(index_beg_spike,1))
        
        temp = intraCell_data (index_beg_spike(ii,1):index_beg_spike(ii,1)+ 80)';
        peak = find(intraCell_data (index_beg_spike(ii,1):index_beg_spike(ii,1)+ 80) == max(temp)) + index_beg_spike(ii,1);
        extracted_spikes_Intra(:,ii+2,i) = intraCell_data(peak-WOI(1)* Ts:peak+WOI(2)* Ts)';

        temp = extraCell_data (index_beg_spike(ii,1):index_beg_spike(ii,1)+ 80)';
        trough = find(extraCell_data (index_beg_spike(ii,1):index_beg_spike(ii,1)+ 80) == min(temp)) + index_beg_spike(ii,1);
        extracted_spikes_Extra(:,ii,i) = extraCell_data (trough-WOI(3)* Ts:trough+WOI(4)* Ts)';
    end
    Header = [file,'ICP',
    %Header
    %Time COl   `


    %Adjust baseline
    temp = -60 - mean(extracted_spikes_Intra(1:10,:,i));
    Baseline = ones(size(extracted_spikes_Intra,1),1)*temp;
    extracted_spikes_Intra(:,:,i) = extracted_spikes_Intra(:,:,i) + Baseline ;
    
    %Scale Amplitude to 1
    temp = 1./max(extracted_spikes_Intra(:,:,i));
    Scalar =  ones(size(extracted_spikes_Intra,1),1)*temp;
    extracted_spikes_Intra(:,:,i) = extracted_spikes_Intra(:,:,i).*Scalar ;
%}
end
figure(3)
plot(extracted_spikes_Extra(:,:,1));
hold on;
%plot(extracted_spikes_Extra(:,:,2));

figure(4)
plot(extracted_spikes_Intra(:,:,1));
hold on;
% add histogram

[PC, SCORE, LATENT, TSQUARE] = princomp(extracted_spikes_Intra');
figure(50)
plotmatrix(SCORE(:,1:5))
[PC, SCORE, LATENT, TSQUARE] = princomp(extracted_spikes_Extra');
figure(51)
plotmatrix(SCORE(:,1:5))

cluster1 = [];
cluster2 = [];

clusters = kmeans(extracted_spikes_Intra',2);

cluster1_ind = find(clusters == 1);
cluster2_ind = find(clusters == 2);

cluster1 = extracted_spikes_Intra(:, cluster1_ind);
cluster2 = extracted_spikes_Intra(:, cluster2_ind);

%
figure(61)

subplot(2,1,1)
plot([1:1:length(cluster1(:,1))]./Ts,cluster1)
ylabel('[mV]')

subplot(2,1,2)
plot([1:1:length(cluster2(:,1))]./Ts,cluster2)
ylabel('[mV]')
xlabel('Time [s]')

figure(62)
plotmatrix(SCORE(cluster1_ind,1:5))
figure(63)
plotmatrix(SCORE(cluster2_ind,1:5))
luster1 = [];
cluster2 = [];

clusters = kmeans(extracted_spikes_Intra',2);

cluster1_ind = find(clusters == 1);
cluster2_ind = find(clusters == 2);

cluster1 = extracted_spikes_Extra(:, cluster1_ind);
cluster2 = extracted_spikes_Extra(:, cluster2_ind);

%
figure(71)

subplot(2,1,1)
plot([1:1:length(cluster1(:,1))]./Ts,cluster1)
ylabel('[mV]')

subplot(2,1,2)
plot([1:1:length(cluster2(:,1))]./Ts,cluster2)
ylabel('[mV]')
xlabel('Time [s]')

figure(72)
plotmatrix(SCORE(cluster1_ind,1:5))
figure(73)
plotmatrix(SCORE(cluster2_ind,1:5))
% %
%plot(extracted_spikes_Intra(:,:,2));



%bas
%extracted_spikes = cat(2,extracted_spikes_Intra,extracted_spikes_Extra);
%   surf(extracted_spikes_Intra(:,1,1),extracted_spikes_Intra(:,2:end,1));
%{
figure(2)
plot(extracted_spikes_Extra,);

extracted_spikes = cat(2,extracted_spikes_Intra,extracted_spikes_Extra);
clear extracted_spikes_Intra;
clear extracted_spikes_Extra;
clear intraCell_data;
clear extraCell_data;
clear data_GT0_shift;
clear data_GT0;

%Adjust baseline
%scale



%Database Structure
%Waveform col
%20 rows header
% Filename
% Protocol (Pulse Burst)
% Channel
% Sweep
% Time in Sweep
% Ts
% Interspike interval
% Baseline Adjustment
% Scale factor

my_beg =  0.095; %0.5015;
my_end =  0.110; %0.530;

chn = [2:N]; % Note: Chn 1 is time

figure(1)
plot(my_data(my_beg*Ts+1:my_end*Ts,1),my_data(my_beg*Ts+1:my_end*Ts,chn))
ylabel('mV')
ylim([-0.02 0.02])





% -------------------------------------------------------------------------



chn = [2:2]; % Note: Chn 1 is time

my_data_proc = my_data;
%artefact = mean(my_data(:,[2:N]),2);

figure(2)
plot(my_data(my_beg*Ts+1:my_end*Ts,1),artefact(my_beg*Ts+1:my_end*Ts))

spike_sort = 1


if spike_sort == 1

figure(4)
plot(my_spikes)
figure(5)
hist(my_ind)


[PC, SCORE, LATENT, TSQUARE] = princomp(my_spikes');
figure(50)
plotmatrix(SCORE(:,1:5))
cluster1 = [];
cluster2 = [];

clusters = kmeans(my_spikes',2)

cluster1_ind = find(clusters == 1);
cluster2_ind = find(clusters == 2);

cluster1 = my_spikes(:, cluster1_ind);
cluster2 = my_spikes(:, cluster2_ind);

%
figure(51)

subplot(2,1,1)
plot([1:1:length(cluster1(:,1))]./Ts,cluster1)
ylabel('[mV]')

subplot(2,1,2)
plot([1:1:length(cluster2(:,1))]./Ts,cluster2)
ylabel('[mV]')
xlabel('Time [s]')

figure(52)
plotmatrix(SCORE(cluster1_ind,1:5))
figure(53)
plotmatrix(SCORE(cluster2_ind,1:5))

% %
% % % -------------------------------------------------------------------------
% %
% %
% %
end

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