% Multielectrode waveform extraction
%   For each electrode threshold and extract spiketimes.
%   Collect waveform from all other electrodes (allow specifying how many
%   and how many samples)

%   Plot the raw data.
%   kmeans cluster
colororder = ['-.r';'-.g';'-.b';'-.c';'-.m';'-.y'];
colororder2 =['r-';'g-';'b-';'c-';'m-';'y-'];
writedirheader = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Data Analysis\Expts\';
% readdirheader =  'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\';
% [readfileheader readfilenumber]=  readaxonDialog(bSliceCell)
% processedFiles = createaxonfilename(readfileheader,readfilenumber(1))
tic
% output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\subset2006_03_20_0009.abf',-1,1);
output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2005_12_13_0004s.abf',-1,1);
toc
% ADD REPORT extracted
figure(1)
plot(output.data(:,2:output.Nchan:end));


%% ADD Preprocessing Filtering


%% CREATE DIR for data
writedirheader = [writedirheader output.sfilename(1:end-4) '\'];
if ~isdir(writedirheader)
    mkdir(writedirheader);
end;

figheader = '';
%%% EXTRACTED from DATA
Ts = 1/output.dt; % Sampling Frequency (Typically 20kHz)
N_chn = output.Nchan; % Number of ADC Channels recorded from
N_samples = output.Nsamples; % Number of ADC Channels recorded from
sw = max(output.xSweeps); % sweeps to set max time

%%% USER INPUT Data parameters
nIntracellsignals =1; % #ADC channels which were Intracellular (Code assumes they are contiguous from 1 to N)
indexoffset = 1  ; % #Offset to ADC channels cause data starts with time column
neighbors = [5;5;5]; % ADC number electrodes that are used with other electrode in spikesorting (waveforms will be concatinated)
% neighbors = -1;

%%% Spike Extraction Parameters
WOI = [750e-6 750e-6; 1500e-6 1500e-6];%Window of interest around Beg_spike -1000us + 2000us
std_threshold = 4; %number of std for threshold
b_intra = 0; % trigger on negative going spike
maxpeakamp = 10000*std_threshold;  %in std
maxpeakwidth =  2e-3;

%%Data output Parameters
breport = 1; % toggle report output
bsave = 1;  % toggle saving data e.g Xspikes,clusters..

index_beg_spike = 0; temp1 = []; temp2 = []; temp2 = [];
% for triggerEE=nIntracellsignals+1:N_chn  % all Extracellular signals , threshold
for triggerEE=4  % all Extracellular signals , threshold
    if ~isinteger(output.data)
        fstd(triggerEE) =  mean(std(output.data(:,indexoffset+triggerEE:N_chn:end)));
    else
        fstd(triggerEE) =  mean(int16std(output.data(:,indexoffset+triggerEE:N_chn:end)));
    end

    fmean(triggerEE) =  mean(mean(output.data(:,indexoffset+triggerEE:N_chn:end)));

    %%% CHANGE code so spikes are in std units
    tic
    %% REmember threshold is mean - std*std_threshold
    temp = int32(FindSpikes2(output.data(:,indexoffset+triggerEE:N_chn:end), (fmean(triggerEE)-(std_threshold*fstd(triggerEE))) ,Ts,b_intra));
    toc
    nSpikes(triggerEE) = size(temp,1);
    index_beg_spike = [index_beg_spike; temp];
    stemp = sprintf('N: %d Spike Detection Parameters\tSTD:\t%1.3f, M:\t%1.3f,\tTH:%d',nSpikes(triggerEE), fstd(triggerEE),fmean(triggerEE),std_threshold)
end
clear temp;
if isempty(find(neighbors == -1))
    for i=1:size(neighbors,2)  %Extract neighbors (column number - offset in data) for triggerEE using ADC numbers
        temp(i)  = find(output.nADCSamplingSeq == neighbors(triggerEE-nIntracellsignals,i));
    end
else
    temp = -1;    %Extract spikes
    %  it doesn't matter which triggerEE is used since spikes are from all
    %  electrodes
    index_beg_spike = unique(index_beg_spike); %% Because there may be spikes detected by multiple electrodes
    % note there is a potential problem because spikes within > 1 sample
    % are allowed
end

[allElectrodesXspikes spikeData temp skipped]  = xSpikesEE2(index_beg_spike,output.data,triggerEE,Ts,N_chn,maxpeakamp+fmean(triggerEE),maxpeakwidth,temp,nIntracellsignals,WOI,fstd);
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
    temp = sprintf('%s%sEEall',spath,figdesc);
    XspikeParams = struct('std',fstd,'mean',fmean,'std_threshold',std_threshold,'maxpeakamp',maxpeakamp,'maxpeakwidth',maxpeakwidth,'WOI',WOI);
    save(temp,'XspikeParams','neighbors','allElectrodesXspikes','spikeData','Spikeskipped','nSpikes');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  xtract INTRAspikes
Intrastd(nIntracellsignals) =  mean(std(single(output.data(:,indexoffset+nIntracellsignals:N_chn:end))));
Intramean(nIntracellsignals) =  mean(mean(single(output.data(:,indexoffset+nIntracellsignals:N_chn:end))));
temp = int32(FindSpikes(output.data(:,indexoffset+nIntracellsignals:N_chn:end), Intramean(nIntracellsignals)+5*Intrastd(nIntracellsignals) ,Ts,1));
INneigh = 0;
[InXspikes InspikeData temp Inskipped]  = xSpikesEE2(temp,output.data,1,Ts,N_chn,maxpeakamp,maxpeakwidth,INneigh,nIntracellsignals,WOI,Intrastd);
clear Inskipped; clear index_beg_spike;
%%%% PLOT DATA
binsize = 2e-3;
range=0.1;
% nbins =        (N_samples/N_chn * sw)/ Ts/(binsize);
nSperBin =  1/output.dt / 1000 * binsize * 1000; % samples per bin
xcorrBinRange = range/binsize;


%%%%%%%%%%%%%%%%
%%  assuming repeated stimulus in Intra ....
%%          Eliminate spikes that occur more then 2 std from mean
%%          then Pick on the first spike from each sweep
temp = std(InspikeData(:,3)); % stdSpikeTimeInSweep
temp1 = mean(InspikeData(:,3));
temp1 = find(InspikeData(:,3)< (temp1+1.5*temp) & InspikeData(:,3)> (temp1-1.5*temp) ); % all spikes less than 2stds from mean spiketime in sweep
i = 0;
for n = 1:max(InspikeData(:,2)) % for all sweeps
    temp2 = min(find(InspikeData(temp1,2)==n)); % find first spike in sweep
    if ~isempty(temp2)
         i = i +1;
         temp(i) = temp2;
    end
end
first_InspikeData = InspikeData(temp,:);

figure(5)
set(gcf,'Position',[500 0 500 1120])
hold off;
subplot(1,3,1)
plot(InXspikes,'k')

axis tight
if ~isempty(InspikeData)
    cross_corr = spiketime_xcorr(InspikeData(:,1),InspikeData(:,1),nSperBin,xcorrBinRange);
    subplot(1,3,2)
    plot(([1:size(cross_corr,2)]-range/binsize)*binsize*1000-1,cross_corr,'-b');% X axis time in Sec
    axis tight
    xlabel('ms');

    subplot(1,3,3) %% plot only first spike
    plot(InXspikes(:,temp),'k');
    axis tight
end

clear InXspikes; clear InspikeData; clear temp;

clear output.data;
% b_intra = 1;
%     %% find Intracellular Spikes
%     temp = int32(FindSpikes(output.data(:,indexoffset+triggerEE:N_chn:end), std_threshold*fstd(triggerEE) ,Ts,b_intra));

% PCA
% if pca_flag ==1
[PC, SCORE, LATENT, TSQUARE] = princomp(single(allElectrodesXspikes'));
% end
% Scatter plot first two PCs
figure(10)
hold off;
set(gcf,'Position',[700 0 700 1120])
subplot(2,2,1)
plot(SCORE(:,1),SCORE(:,2),'.')
axis tight;
title('PCA1 vs PCA2');
subplot(2,2,2)
plot(LATENT,'x');
axis tight;
title('EigenValues');
subplot(2,2,3:4)
imagesc(hist3([SCORE(1:end,1) SCORE(1:end,2)],[200 200]));
title('PCA1 vs PCA2');

if breport
    dirtemp = 'REPORT';
    figdesc = 'PCA12';
    spath = [writedirheader dirtemp '\'];
    if (~isdir(spath))
        mkdir(spath);
    end
    temp = sprintf('%s%sEE%d.png',spath,figdesc,triggerEE);
    %     print('-dtiffn',temp);
    print('-dpng',temp);

end
clear PC; clear LATENT; clear TSQUARE;
%number of principle components to use
nPC = 20;
% Overcluster data
num_c = 16;
tt = size(colormap,1);
temp = colormap;
my_colors=temp([1:round(tt/num_c):tt],:);
my_colors = [my_colors; my_colors];
my_style =['.','o','x','+','*','s','d','v','^','<','>','p','h','.','o','x'
    ,'.','o','x','+','*','s','d','v','^','<','>','p','h','.','o','x'];
my_syle = [my_style, my_style, my_style, my_style];

% kmeans_flag = 1;
% if kmeans_flag == 1
tic
%  [klusters Center]= kmeans(single(allElectrodesXspikes)',num_c);
[klusters Center]= kmeans(SCORE(1:end,1:nPC),num_c);
toc
% end
% Show klusters in 2dim PCA space
figure(20)
hold off
set(gcf,'Position',[700 0 1500 1120]);
subplot(2,1,1)
for i =1:num_c
    pp=plot(SCORE(klusters==i,1),SCORE(klusters==i,2),'.');
    set(pp,'color',my_colors(i,:));
    % set(pp,'Marker',my_style(i));
    hold all
    %     i
    %     pause;
end
legend('1', '2', '3' ,'4' ,'5', '6', '7', '8', '9',...
    '10','11','12','13','14','15','16');
title(['Overclustered in PCA1 vs PCA2 (nPC:' num2str(nPC) ')']);

% Combine klusters
% Y = pdist(Center);
Y = pdist(Center,'mahalanobis');
Z = linkage(Y);  % cluster function might be used to autocluster
% dendrogram(Z)
t = .5*(max(Z(:,3)));
subplot(2,1,2);
dendrogram(Z,0,'colorthreshold',t,'orientation','top'); % to show all klusters
%% ADD number of spikes in each cluster
if bsave
    dirtemp = '';
    figdesc = 'ClusterData';
    spath = [writedirheader dirtemp '\'];
    if (~isdir(spath))
        mkdir(spath);
    end
    temp = sprintf('%s%sEE%d',spath,figdesc,triggerEE);
    %     print('-dtiffn',temp);
    save(temp,'klusters','Center','num_c')  %% NOTE font size too big on axis (overlaps)
end
if breport
    dirtemp = 'REPORT';
    figdesc = 'Cluster';
    spath = [writedirheader dirtemp '\'];
    if (~isdir(spath))
        mkdir(spath);
    end
    temp = sprintf('%s%sEE%d.png',spath,figdesc,triggerEE);
    %     print('-dtiffn',temp);
    print('-dpng',temp)  %% NOTE font size too big on axis (overlaps)
end

%% NEED To add report of data
clear SCORE;
pack;
n_units = 10;
tic
T = cluster(Z,'maxclust',n_units) ;
toc

unit_Nspikes = zeros(1,n_units);
NN = 4; % subplots per figure
NP = 5 % plots per cluster;
Nwindows = ceil(n_units/NN)
for k = 1:1:Nwindows % Number of figures
    figure(k+300)
    set(gcf,'Position',[200 0 (400+NP*200) 1120])
    for tt = 1:NN %% Number of subplots per figure
        j = (k-1)*NN + tt;
        tsubp = (tt-1)*NP+1;
        subplot(NN,NP,tsubp:tsubp+1)
        %        single_unit{j} = unit{j};
        single_unit{j} = int32(find(T==j));
        unit_Nspikes(j) = 0;  temp=[];
        for i=1:length(single_unit{j}) % Sequentially plots spikes from each cluster in a Unit
            if i ==1
                hold off;
            end
            temp = [temp allElectrodesXspikes(:,klusters==single(single_unit{j}(i)))];
            %% Extract spiketimes and sweep
            temp_spikeData{j} = spikeData(klusters == single_unit{j}(i),1:3);

        end
        spikes_to_plot{j} = temp;
        unit_Nspikes(j) = unit_Nspikes(j) + size(temp,2);

        plot([1:1:size(allElectrodesXspikes,1)].*Ts,spikes_to_plot{j},'k')
        axis tight
        hold on
        m= min([size(single_unit{j},1) 5]);
        temp = sprintf('Unit:%d Spikes:%d C:%s',j,unit_Nspikes(j),num2str(single_unit{j}(1:m),'%d,'));
        title(temp)


        subplot(NN,NP,tsubp+2)
        binsize = 1e-3; %(sec)
        range = .04;% range(sec) of to compute cross correlation

        title([num2str(temp) ' spikes in cluster'])
        xlabel('[sec]')
        ylabel('[mV]')
        unit_spiketimesX = single(spikeData(klusters==j,1));

        nSperBin =       1/ output.dt/1000 * binsize * 1000;
        xcorrBinRange = range/binsize;
        cross_corr = spiketime_xcorr(unit_spiketimesX,unit_spiketimesX,nSperBin,xcorrBinRange);
        plot(([1:size(cross_corr,2)]-range/binsize)*binsize*1000-1,cross_corr,'-g');% X axis time in Sec
        axis tight
        ylim([0 .1])

        %% ADD spike rate, EI?,

        %% PLOT average Intracellular signal
        subplot(NN,NP,tsubp+3)
        if nIntracellsignals >0
            % correct index
            for jjj=1:nIntracellsignals
                temp = output.data(:,1+jjj:output.Nchan:end);
                temp1 = zeros(size(temp_spikeData{j},1),int32(Ts*50e-3)+1);
                for jjjj=1:size(temp_spikeData{j},1)
                    temp2 = (temp_spikeData{j}(jjjj,1));
                    if temp2
                    temp1(jjjj,:) = temp(temp2:(temp2 + int32(Ts*50e-3)));
                    end
                end
                temp_Intra{j} = temp1;
                hold on;
                plot([1:size(temp1(jjjj,:),2)]*1e3/Ts,single(temp1)'*output.gain(jjj),'-k');
                %                        plot([1:size(temp1(jjjj,:),2)]*1e3/Ts,(mean(temp1,1)+std(single(temp1),1))*output.gain(jjj),'--k');
                plot([1:size(temp1(jjjj,:),2)]*1e3/Ts,mean(temp1,1)*output.gain(jjj),colororder2(jjj,:));
                axis tight
                xlabel('ms')
            end
        end
        
        subplot(NN,NP,tsubp+4)
        binsize = 5e-3;
        range=0.1;
%         nbins =        (N_samples/N_chn * sw)/ Ts/(binsize);
        nSperBin =       1/ output.dt/1000 * binsize * 1000;
        xcorrBinRange = range/binsize;
        unit_spiketimesX = single(spikeData(klusters==j,1));
        cross_corr = spiketime_xcorr(first_InspikeData(:,1),unit_spiketimesX,nSperBin,xcorrBinRange);
        plot(([1:size(cross_corr,2)]-range/binsize)*binsize*1000-1,cross_corr,'-b');% X axis time in Sec
        axis tight


    end


end
%
%%% ACCEPT/JOIN
accept_units{1} = {1}; %% this number is the cluster from Cluster (NOT kmeans)
accept_units{2} = {3};
accept_units{3} = {5,14};
accept_units{4} = {6};
accept_units{5} = {10};
accept_units{6} = {11};
accept_units{7} = {12};
accept_units{8} = {15};
accept_units{9} = {9};
%% REJECT
reject_unit = [16:20 7 8];

%% Extract Units
good_units = struct([]);
for j = 1:size(accept_units,2)
    clear temp; clear temp1;
    for jj = 1:size(accept_units{j},2)
        if jj ==1
            temp = spikes_to_plot{cell2mat(accept_units{j}(1,jj))};
            temp1 = temp_spikeData{cell2mat(accept_units{j}(1,jj))};
            temp2 = single_unit{cell2mat(accept_units{j}(1,jj))};
        else
            temp = [temp spikes_to_plot{cell2mat(accept_units{j}(1,jj))}];
            temp1 = [temp1; temp_spikeData{cell2mat(accept_units{j}(1,jj))}];
            temp2 = [temp2; single_unit{cell2mat(accept_units{j}(1,jj))}];

        end

    end
    good_unit(j).spikes = temp;
    good_unit(j).spikeData = temp1;
    if j ==1
        removeKluster = temp2;
    else
        removeKluster = [removeKluster; temp2];
    end
end
removeKluster = [removeKluster; reject_unit']';

removeSpike = [];keepSpike = [];
for i = 1: size(klusters,1)
    if any(removeKluster == klusters(i))
        removeSpike=[removeSpike i];
    else
        keepSpike = [keepSpike i];
    end
end

% recluster w/o
[PC, SCORE, LATENT, TSQUARE] = princomp(single(allElectrodesXspikes(:,keepSpike)'));





%MANUAL
unit = cell(n_units,1);
unit{1} = [5 23 32]; % A
unit{1} = [10 6 24 60 7] % A
unit{2} = [19 27 63 55 62]%B
unit{3} = [15 4 2 28]%C
unit{4} = [1]%D
unit{5} = [3]%E
unit{6} = [6]%F
unit{7} = [7]%G
unit{8} = [8]%H
unit{9} = [9]%I
unit{10} = [10]%J
unit{11} = [12]%K
unit{12} = [13]%L
unit{13} = [14]%M
unit{14} = [39 52 54 48 44 50 29 47 41 22 25 36 51 53 5 21]%N
% % unit{16} = []%O
% % unit{15} = [58 62]%P
% %
%

basis_hist = [0.110:0.0005:0.150]%save electrode
%%%%%%%%%%%%%%%%%%%% CODE for selecting units
%PLOT cells together (display mean and std)
flagged_units = [4 7];     % MANUAL ENTRY
figure(200);
hold off
for i = 1:size(selected_units,2)
    single_unit = find(T==flagged_units(i));
    temp = mean(allElectrodesXspikes(:,klusters==single_unit)');
    plot(temp,colororder2(i,:));
    hold on;
    plot(temp-std(allElectrodesXspikes(:,klusters==single_unit)'),colororder(i,:));
    plot(temp+std(allElectrodesXspikes(:,klusters==single_unit)'),colororder(i,:));
end


%%CROSS-CORRELATION
selected_units = [3 4 6 5 7 9];
binsize = 1e-3; %(sec)
range = .1;% range(sec) of to compute cross correlation
hold off
for k = 1:size(selected_units,2) % Number of figures (1 for each cell)
    figure(k+500);
    hold off
    set(gcf,'Position',[700 0 700 1120])

    for j = 1:size(selected_units,2) %% Number of subplots per figure (1 for each cell)
        subplot(size(selected_units,2),1,j)
        if j == 1;
            hold off;
        end

        unit_spiketimesX = single(spikeData(klusters==selected_units(k),1));
        unit_spiketimesY = single(spikeData(klusters==selected_units(j),1));
        % Bin Data
        %         binned_spikesX = bin(unit_spiketimesX(:,1),(N_samples/N_chn * sw)/ Ts/(binsize))/(binsize);
        %         binned_spikesY = bin(unit_spiketimesY(:,1),(N_samples/N_chn * sw)/ Ts/(binsize))/(binsize);
        %         %         plot((1:size(binned_spikes,2))*(binsize),binned_spikes(1,:))  % X axis time in Se
        %
        %         cross_corr = xcorr(binned_spikesX,binned_spikesY,range/binsize,'coeff');
        nSperBin =       1/ output.dt/1000 * binsize * 1000;
%         nbins =        (N_samples/N_chn * sw)/ Ts/(binsize);
        xcorrBinRange = range/binsize;
        cross_corr = spiketime_xcorr(unit_spiketimesX,unit_spiketimesY,nSperBin,xcorrBinRange);
        %% may not be centered right
        plot(([1:size(cross_corr,2)]-range/binsize)*binsize*1000-1,cross_corr,'ob');% X axis time in Sec
        axis tight
        xlabel('ms');
        %         temp = sprintf('%s%d-%d (%d,%d) Bin:%1.2gms',figheader,k,j,unit_Nspikes(k),unit_Nspikes(j),binsize*1000);
        temp = sprintf('%s%d-%d Spikes:%d,%d, Bin:%1.2gms',figheader,k,j,1,1,binsize*1000);
        title(temp);
    end
    %% save data image to disk
    breport = 0;
    if (breport)
        dirtemp = 'REPORT';
        figdesc = 'XCorr';
        spath = [writedirheader dirtemp '\'];
        if (~isdir(spath))
            mkdir(spath);
        end
        temp = sprintf('%s%s_Unit_%d.png',spath,figdesc,k);
        %     print('-dtiffn',temp);
        print('-dpng',temp);
    end
end
%%%%%%%%%% SELECT UNITS
% selected_units = [3 4 7 8];
%     binsize = 1e-3; %(sec)
%     range = .1;% range(sec) of to compute cross correlation
%%AUTOCORRELATION
% figure(300);
% hold off
% sw = output.xSweeps; % sweeps to set max time
%  for j = 1: size(selected_units,2)
%     unit_spiketimes = spikeData(klusters==selected_units(j),1);
% % Bin Data
%     binned_spikes = bin(unit_spiketimes(:,1),(N_samples/N_chn * sw)/ Ts/(binsize))/(binsize);
%     plot((1:size(binned_spikes,2))*(binsize),binned_spikes(1,:))  % X axis time in Sec
%     hold on
%     cross_corr = xcorr(binned_spikes,range/binsize,'coeff');
%     if j >1
%         pause;
%     end
%     plot(([1:size(cross_corr,2)]-range/binsize)*binsize,cross_corr,'-b');
%     axis tight
%
% %     temp = sprintf('Unit:%d Spikes:%d',j,unit_Nspikes(j));
%     title([num2str(selected_units(j)) 'AutoCorrelation'])
% end
%%%%%%%%%%%%%%%%%%%%%%%
%
% my_times = cell(n_units,1);
% my_sweeps = cell(n_units,1);
% NN = 4;
% Nwindows = ceil(n_units/NN)
% for k = 1:1:Nwindows % Number of figures
%    figure(k+100)
%    hold off
%    set(gcf,'Position',[420 34 816 910])
%    for tt = 1:NN %% Number of subplots per figure
%        j = (k-1)*NN + tt;
%        subplot(NN,1,tt)
%        single_unit = unit{j};
%        single_unit = find(T==j);
%
%        temp = 0;
%        for i=1:length(single_unit) % Sequentially plots spikes from each cluster in a Unit
%            unit_spike_ind(i) = klusters==single_unit(i);
%            spikes_to_plot = allElectrodesXspikes(:,klusters==single_unit(i));
% %            my_times{j} = [my_times{j}
% % spike_time(klusters==single_unit(i))] ;
% %            my_sweeps{j} = [my_sweeps{j}
% % sw_index(klusters==single_unit(i))]
% %            temp = temp + size(spikes_to_plot,2);
%
%            plot([1:1:size(spikes_to_plot,1)].*Ts,spikes_to_plot,'k')
%            axis tight
%            % ylim([-0.07 0.05])
%            axis off
%            title(num2str(j))
%            %text(1e-4,0,num2str(size(spikes_to_plot,2)))
%            hold on
%        end
%        % title([num2str(temp) ' spikes in cluster'])
%        % xlabel('[sec]')
%        % ylabel('[mV]')
%    end
%
%
% end
% basis_hist = [0.110:0.0005:0.150]%save electrode
%
%
% %PLOT cells together
% selected_units = [7 8];
% for i = 1:size(selected_units,2)
%     figure(200);
%     plot(allElectrodesXspikes(:,selected_units(i)));
%        set(pp,'color',my_colors(i,:));
%        hold on;
% end
%