% Multielectrode waveform extraction
%   For each electrode threshold and extract spiketimes.
%   Collect waveform from all other electrodes (allow specifying how many
%   and how many samples)

%   Plot the raw data.
%   kmeans cluster
writedirheader = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Data Analysis\Expts\';
% readdirheader =  'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\';
% [readfileheader readfilenumber]=  readaxonDialog(bSliceCell)
% processedFiles = createaxonfilename(readfileheader,readfilenumber(1)) 
tic
output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_03_20_0008.abf',-1,1);
toc
% REPORT extracted

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
%neighbors = [8 11;5 11; 5 8]; % ADC number electrodes that are used with other electrode in spikesorting (waveforms will be concatinated)
neighbors = -1;

%%% Spike Extraction Parameters
WOI = [1000e-6 1000e-6; 1500e-6 1500e-6];%Window of interest around Beg_spike -1000us + 2000us
std_threshold = -5; %number of std for threshold
b_intra = 0; % trigger on negative going spike
maxpeakamp = 10000*std_threshold;  %in std
maxpeakwidth =  2e-3;  

%%Data output Parameters
breport = 1; % toggle report output
bsave = 1;  % toggle saving data e.g Xspikes,clusters..

index_beg_spike = 0;
for triggerEE=nIntracellsignals+1:N_chn  % all Extracellular signals , threshold
    fstd(triggerEE) =  mean(std(single(output.data(:,indexoffset+triggerEE:N_chn:end))));
    fmean(triggerEE) =  mean(mean(single(output.data(:,indexoffset+triggerEE:N_chn:end))));

    temp = int32(FindSpikes(output.data(:,indexoffset+triggerEE:N_chn:end), std_threshold*fstd(triggerEE) ,Ts,b_intra));
    nSpikes(triggerEE) = size(temp,1);
    index_beg_spike = [index_beg_spike; temp];
     stemp = sprintf('N: %d Spike Detection Parameters\tSTD:\t%1.3f, M:\t%1.3f,\tTH:%d',nSpikes(triggerEE), fstd(triggerEE),fmean(triggerEE),std_threshold)
end
    index_beg_spike = unique(index_beg_spike);
    %Extract spikes
    %  it doesn't matter which triggerEE is used since spikes are from all
    %  electrodes
    [allElectrodesXspikes spikeData temp skipped]  = xSpikesEE2(index_beg_spike,output.data,triggerEE,Ts,N_chn,maxpeakamp,maxpeakwidth,neighbors,nIntracellsignals,WOI,fstd);
    for jj = 1:4
        Spikeskipped(jj) = sum(skipped == -1*jj);
    end
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
% b_intra = 1;
%     %% find Intracellular Spikes
%     temp = int32(FindSpikes(output.data(:,indexoffset+triggerEE:N_chn:end), std_threshold*fstd(triggerEE) ,Ts,b_intra));

% PCA
% if pca_flag ==1
[PC, SCORE, LATENT, TSQUARE] = princomp(single(allElectrodesXspikes'));
% end
plot(LATENT,'x');
% Scatter plot first two PCs
figure(10)
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
imagesc(hist3([SCORE(1:end,1) SCORE(1:end,2)],[400 400]));
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

%number of principle components to use
nPC = 30;
% Overcluster data
num_c = 64;
tt = size(colormap,1);
temp = colormap;
my_colors=temp([1:round(tt/num_c):tt],:);
my_colors = [my_colors; my_colors];

% kmeans_flag = 1;
% if kmeans_flag == 1
tic
% [clusters Center]= kmeans(single(allElectrodesXspikes)',num_c);
[clusters Center]= kmeans(SCORE(1:end,1:nPC),num_c);
toc
% end
% Show clusters in 2dim PCA space
figure(20)
hold off
set(gcf,'Position',[1500 0 1500 1120]);
subplot(2,1,1)
for i =1:num_c
    pp=plot(SCORE(clusters==i,1),SCORE(clusters==i,2),'.');
    set(pp,'color',my_colors(i,:));
    % set(pp,'Marker',my_style(i));
    hold all
end
legend('1', '2', '3' ,'4' ,'5', '6', '7', '8', '9',...
    '10','11','12','13','14','15','16');
title(['Overclustered in PCA1 vs PCA2 (nPC:' num2str(nPC) ')']);

% Combine clusters
% Y = pdist(Center);
Y = pdist(Center,'mahalanobis');
Z = linkage(Y);  % cluster function might be used to autocluster 
% dendrogram(Z)
 t = .5*(max(Z(:,3))); 
subplot(2,1,2);
dendrogram(Z,0,'colorthreshold',t,'orientation','top'); % to show all clusters
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
    save(temp,'clusters','Center','num_c')  %% NOTE font size too big on axis (overlaps)
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

n_units = 20;
tic
T = cluster(Z,'maxclust',n_units) ;
toc
unit = cell(n_units,1);
unit{1} = [18 49 59 35] % A
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
my_times = cell(n_units,1);
my_sweeps = cell(n_units,1);
unit_Nspikes = zeros(1,n_units);

NN = 4; % subplots per figure
Nwindows = ceil(n_units/NN)
    binsize = 1e-3; %(sec)
    range = .01;% range(sec) of to compute cross correlation
for k = 1:1:Nwindows % Number of figures
   figure(k+200)
   set(gcf,'Position',[700 0 700 1120])
   for tt = 1:NN %% Number of subplots per figure
       j = (k-1)*NN + tt;
       tsubp = (tt-1)*3+1;
       subplot(NN,3,tsubp:tsubp+1)
       single_unit = unit{j};
%        single_unit = find(T==j);
        unit_Nspikes(j) = 0;
       for i=1:length(single_unit) % Sequentially plots spikes from each cluster in a Unit
           if i == 1;
               hold off;
           end
           spikes_to_plot = allElectrodesXspikes(:,clusters==single_unit(i));
%            my_times{j} = [my_times{j}
% spike_time(clusters==single_unit(i))] ;
%            my_sweeps{j} = [my_sweeps{j}
% sw_index(clusters==single_unit(i))]
           unit_Nspikes(j) = unit_Nspikes(j) + size(spikes_to_plot,2);

           plot([1:1:size(spikes_to_plot,1)].*Ts,spikes_to_plot,'k')
           axis tight
%            axis off
           %text(1e-4,0,num2str(size(spikes_to_plot,2)))
           hold on
       end
       m= min([size(single_unit,2) 5])
       temp = sprintf('Unit:%d Spikes:%d C:%s',j,unit_Nspikes(j),num2str(single_unit(1:m)));
       title(temp)
       
       subplot(NN,3,tsubp+2)
       % title([num2str(temp) ' spikes in cluster'])
       % xlabel('[sec]')
       % ylabel('[mV]')
%         unit_spiketimesX = single(spikeData(clusters==j,1));
% 
%                nbins =        (N_samples/N_chn * sw)/ Ts/(binsize);
%         xcorrBinRange = range/binsize;
%         cross_corr = spiketime_xcorr(unit_spiketimesX,unit_spiketimesX,nbins,xcorrBinRange);
%         plot(([1:size(cross_corr,2)]-range/binsize)*binsize*1000-1,cross_corr,'.b');% X axis time in Sec
%         axis tight
%         ylim([0 .1])
        
        %% ADD spike rate, EI?, 
        
   end


end

basis_hist = [0.110:0.0005:0.150]%save electrode 
%%%%%%%%%%%%%%%%%%%% CODE for selecting units
%PLOT cells together (display mean and std)
flagged_units = [4 7];     % MANUAL ENTRY
colororder = ['-.r';'-.g';'-.b';'-.c';'-.m';'-.y'];
colororder2 =['r-';'g-';'b-';'c-';'m-';'y-'];
figure(200);
hold off
for i = 1:size(selected_units,2)
    single_unit = find(T==flagged_units(i));
    temp = mean(allElectrodesXspikes(:,clusters==single_unit)');
    plot(temp,colororder2(i,:));
    hold on;
    plot(temp-std(allElectrodesXspikes(:,clusters==single_unit)'),colororder(i,:));
    plot(temp+std(allElectrodesXspikes(:,clusters==single_unit)'),colororder(i,:));
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
 
        unit_spiketimesX = single(spikeData(clusters==selected_units(k),1));
        unit_spiketimesY = single(spikeData(clusters==selected_units(j),1));
        % Bin Data
%         binned_spikesX = bin(unit_spiketimesX(:,1),(N_samples/N_chn * sw)/ Ts/(binsize))/(binsize);
%         binned_spikesY = bin(unit_spiketimesY(:,1),(N_samples/N_chn * sw)/ Ts/(binsize))/(binsize);
%         %         plot((1:size(binned_spikes,2))*(binsize),binned_spikes(1,:))  % X axis time in Se
%         
%         cross_corr = xcorr(binned_spikesX,binned_spikesY,range/binsize,'coeff');
        nbins =        (N_samples/N_chn * sw)/ Ts/(binsize);
        xcorrBinRange = range/binsize;
        cross_corr = spiketime_xcorr(unit_spiketimesX,unit_spiketimesY,nbins,xcorrBinRange);
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
%     unit_spiketimes = spikeData(clusters==selected_units(j),1); 
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
%            unit_spike_ind(i) = clusters==single_unit(i);
%            spikes_to_plot = allElectrodesXspikes(:,clusters==single_unit(i));
% %            my_times{j} = [my_times{j}
% % spike_time(clusters==single_unit(i))] ;
% %            my_sweeps{j} = [my_sweeps{j}
% % sw_index(clusters==single_unit(i))]
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