function [usf unitdatafile spikes] =  sortSpikes(datafilename,NsortAlg,savefilename,bUSERMANUAL,bCrossCorr)
% function [usf unitdatafile spikes] =  sortSpikes(datafilename,NsortAlg,savefilename,bUSERMANUAL,bCrossCorr)
% NsortAlg = 2, KlustaKwik
% NsortAlg = 1, Chronux Kmeans 
% bUSERMANUAL = 1 (only in with chronux sorting) user can manually define
% which kluster form a unit
%
% OUTPUT:
%      usf - savefilename for unit data
%      unitdatafile - struct with data
%      spikes - struct of chronux klustering (is empty for klustakwik
%      cluster)
%% Load in spiketimes file
% (generated with extractDAQ_Spikes.m)
%now you have several options:
% KMEANS: s_kmeanssort.m or Chronux SpikeSortingDemo.m
% KlustaKwik (EM)
% ICA first then KMEANS
% the code to run each of these is below
%
% so far (with one data set,  Kmeans on entire waveform then manually joined seems to looks the
% best in terms of agreeing with what I would do by eye.
% BA
if nargin < 4
bUSERMANUAL = 0; % User manually defines klusters
end
if ~exist('bCrossCorr')
bCrossCorr = 0; % cross correlation spike times of units
end

bICA = 1;
load(datafilename)

% parameters
bsave = 1; % save unit data to file
% prepare output struct
swlength = datafile(1).maxsample;
dt = datafile.dt;
%% concatinate waveforms,spiketimes
[stimes sw swstimes conWV meanWV] = concSpikeData(datafile);

%             [pkpk] = getpkpk(datafile,0); % USE to sort by amplitude on
%             each electrodes
%% plot range and amplitude
figure(30);clf;
N = length(datafile(1).spiketimes);
for i = 1: N
    st{i} = num2str(i);
    
    AMP = min(conWV((end/N)*(i-1)+1:(end/N)*(i),:));
    [a x] = hist(AMP,100);
    %     subplot(N,1,i)
    subplot(2,1,1);
    stairs(x,a); hold all;
    title('Trough amp')
    
    ran = range(conWV((end/N)*(i-1)+1:(end/N)*(i),:));
    [a x] = hist(ran,100);
    subplot(2,1,2);
    stairs(x,a); hold all;
    xlabel('mV')
    title('range')    
end

legend(st);

    
%%
% thresholds% *** ASSUME that threshold doesn't change across files 
for i = 1: length(datafile(1).spiketimes); unitdatafile.th(i) = datafile(1).spiketimes(i).th; end
 
unitdatafile.filename = {datafile.filename};
unitdatafile.MAXTRIGGER = {datafile.MAXTRIGGER};
unitdatafile.maxsample = {datafile.maxsample};
unitdatafile.params = {datafile.params};
unitdatafile.LOADPATH = datafile(1).LOADPATH;
clear datafile;

%% get filename and path info
indtemp = strfind(datafilename,'\');FILENAME = datafilename(indtemp(end)+1:end);
PATH = datafilename(1:indtemp(end));

sAnn = FILENAME;
%% deal with outliers using CHX
clear spikes
% spikes.waveforms = double((conWV-meanWV)');
spikes.waveforms = double((conWV)');

spikes.spiketimes = stimes*dt;
spikes.Fs = 1/dt;
spikes.threshT= 50; %round(WOI(1)/1e3/dt)+1+50;
% if bZSCORE; spikes.threshV = [-TH_STD inf]; else spikes.threshV = [-TH_STD*mean(STDEV) inf]; end
spikes.threshV = [mean(unitdatafile.th) inf];
spikes.thresh = spikes.threshV(1); % this seems to have been added in Dan Hill's version
%***** run Chronux sorting (kmean)
% remove outliers
spikes = ss_outliers(spikes);
figure(3);
subplot(2,1,1);plot(spikes.waveforms'); axis tight; title('Centered Data w/ Outliers Removed');
subplot(2,1,2);plot(spikes.outliers.waveforms'); title('Outliers');axis tight;

switch NsortAlg
    
    case 1
        sortingAlg = 'CHX';
        %********** Chronux sort % see SpikeSortingDemo.m for details
        clear  unit;
        %         bICA = 0; % doesn't seem to work well for this clustering
        %         if bICA;    [icasig A W] = FASTICA (conWV-meanWV, 'numOfIC', 5,'verbose', 'on', 'displayMode', 'on'); spikes.waveforms = icasig';
        %         else spikes.waveforms = double((conWV-meanWV)');end
        
        spikes = ss_kmeans(spikes);
        figure(6);set(gcf,'Position',[ 2282        -104         983        1070]);
        subplot(6,2,[2 4]);
        ssg_databrowse2d(spikes);
        % ssg_databrowse3d(spikes);
        subplot(6,2,[6 12]);clusterXT(spikes, [],[],[],1); title('Final Clusters');
%         subplot(6,2,[6 12]);clusterXT(spikes); title('Final Clusters');
        
        % try and use chronux automatic grouping
        spikes = ss_energy(spikes); spikes = ss_aggregate(spikes);
        
        figure(6);subplot(6,2,[1 3]);
        ssg_databrowse2d(spikes);
        
        
        set(gcf, 'Renderer', 'OpenGL');
        subplot(6,2,[5 7]);clusterXT(spikes, spikes.hierarchy.assigns,[],[],0); title('Final Clusters');
%         subplot(6,2,[5 7]);clusterXT(spikes, spikes.hierarchy.assigns); title('Final Clusters');
        subplot(6,2,9);aggtree(spikes); title('Aggregation Tree');
        
        % add plot s stationarity (convolve with gaussian) and overlay lines
        subplot(6,2,11); plot(spikes.spiketimes, spikes.hierarchy.assigns, '.'); xlabel('Time (sec)'); ylabel('Cluster #');axis tight;
        %% extract "isolated" units
        %         K = unique(spikes.hierarchy.assigns);
        %         K=K(find(K>0));
        %
        %
        %         NK = length(K)
        %         for i=1: % Get number of spikes in each cluster
        %             I = find(spikes.hierarchy.assigns==K(i));
        %             K_N(i) = length(I);
        %         end
        %         PICK_KLUSTERS = find(K_N>50);
        assign = spikes.hierarchy.assigns; % default (may be changed if user specifies cluster by hand below)
        
        if bUSERMANUAL
            sortingAlg = 'MANCHX';
            
            s = 0;i =1;bloop = 1;
            while bloop % USER enters kluster info (1 unit at a time with spaces seperationg cluster number)
                ss = sprintf('Unit%d, Enter Kluster #s (Quit-q): ',i);
                s = input(ss,'s');
                if isempty(s) | ~isempty(strfind(lower(s),'q')); % end
                    bloop  =0; 
                else
                    UserK(i).k = str2num(s);
                    i=i+1;
                end
            end
            
            if i>1 % if user entered any clusters
                oldassign = spikes.overcluster.assigns;
                assign = zeros(size(oldassign)); % note zero is ungrouped
                for i = 1:length(UserK)
                    ind = [];
                    for j  = 1:length(UserK(i).k)
                        ind = [ind; find(oldassign== UserK(i).k(j))];
                    end
                    assign(ind) = i;
                end
            end
                    
         end
        
        figure(8); set(gcf,'Position',[1689         559         560         420]); colormap jet; % THis plot overlaps with plotKluster but has spiking over time so I kept it
        showclust(spikes, spikes.hierarchy.assigns);
% figure(9);
% set(gcf,'Position',[710        -145         560         420])
% separation_analysis(spikes)

        [ifid PICK_KLUSTERS]=local_plotclusters(assign, conWV, stimes,dt,sAnn,bCrossCorr) ;
        clear conWV meanWV
        
        if ~exist('unit'); unitL = 0; else L = length(); end;
        Kluster2extract = PICK_KLUSTERS; % can extract a subset of all clusters
        for i=1:length(Kluster2extract) % extract
            I = find(assign==Kluster2extract(i));
            unit(unitL+i).waveforms = []; % spikes.waveforms(I);
            unit(unitL+i).spiketimes = spikes.spiketimes(I);
            unit(unitL+i).sw = sw(I);
            unit(unitL+i).sw_spiketimes = swstimes(I)*dt;
            unit(unitL+i).index_in_spikes = I;
        end
        
        
        %%% ***************** TO DO
        %             MAKE kluster plotting (below) work with CHX or KKW klusters
        %%%
        
        
        % %%  ** USEFUL remove extracted unit's spikes from data set  (now can rerun Chronux  sorting (above)
        % nI = ones(size(stimes)); for i = 1:length(unit);   nI(unit(i).index_in_spikes) = 0; end;
        % nI = find(nI);
        % clear spikes;
        % spikes.waveforms = double(conWV(:,nI)');
        % spikes.spiketimes = stimes(nI);
        % spikes.Fs = 1/dt;
        % spikes.threshT= 50; %round(WOI(1)/1e3/dt)+1+50;
        
        %%
        
        
    case 2 %% TOD MUST REMOVE OUTLIERS
        
        % remove outliers
        conWV = conWV(:,spikes.outliers.goodinds);
        meanWV = meanWV(:,spikes.outliers.goodinds);
        stimes = stimes(spikes.outliers.goodinds);
        
        spikes = [] ; % this struct is only defined for CHRONUX clustering
        clear K_N unit;
        sortingAlg = 'KKW';
        if bICA
            [icasig A W] = FASTICA (conWV-meanWV, 'numOfIC', 5,'verbose', 'on', 'displayMode', 'on');
            [assign nk] = kkwik(icasig',[PATH '\KKW\' FILENAME]);
        else
            %             [assign nk] = kkwik(pkpk',[PATH '\KKW\' FILENAME]); % to sort
            %             by pkpk amplitude on each electrode
            
            [assign nk] = kkwik((conWV-meanWV)',[PATH '\KKW\' FILENAME]);
        end
        %
        
        
        [ifid PICK_KLUSTERS]=local_plotclusters(assign, conWV, stimes,dt,sAnn,bCrossCorr);
        clear conWV meanWV;
        
        %ADD JOIN clusters...
        
        Kluster2extract = PICK_KLUSTERS; % can extract a subset of all clusters
        for i=1:length(Kluster2extract)
            I = find(assign==Kluster2extract(i));
            unit(i).waveforms = []; %% NOT SAVING Waveform again (it can be loaded from save extracted waveforms)
            unit(i).spiketimes = stimes(I)*dt;
            unit(i).sw = sw(I);
            unit(i).sw_spiketimes = swstimes(I)*dt;
            unit(i).index_in_spikes = I;
        end
 end

clear  swstimes sw;

unitdatafile.dt = dt; %unit(1).swlength = maxsample;
unitdatafile.unit = unit;
unitdatafile.bICA=  bICA;

if exist('savefilename')
    sf = sprintf('%s_SSFig_%s',FILENAME,sortingAlg);
    set(0,'CurrentFigure',ifid+1);    savefigure(PATH,'',sf,'eps');
    
    usf = sprintf('%s%s_SSunitdata_%s',PATH,FILENAME,sortingAlg);save(usf,'unitdatafile');
end
clear  unit;


function [ifid PICK_KLUSTERS]=local_plotclusters(assign, conWV, stimes,dt,sAnn,bCrossCorr)
% plots klusters with greater than X spikes
% ifid (figure number of plots)
% PICK_KLUSTERS klusters that are ploted

ifid = 4;
NK = unique(assign); % plot
NK = NK(find(NK>0)); % 0 is outliers

% don't bother with klusters with less than 50 spikes
clear K_N;
for i= 1:length(NK); K_N(i) = sum(assign == NK(i)); end
PICK_KLUSTERS = NK(find(K_N>20));  % only take cluster with more than X spikes
display('Only take Klusters with > 20 spikes');

if 1 % *** analyze klusters
    figure(ifid); clf;   set(gcf, 'Renderer','painters','Position',[2414         493         867         488]);
    figure(ifid+1); clf;set(gcf, 'Renderer','painters','Position',[ 1683        -125         966        1098]);
    temp = zeros(1,max(stimes),'int16'); % for creating spikerate as a function of time
    MAXSPIKESPLOT = 500; YRANGE = [1000 -1000]; clear  h; % number of spikes in each cluster
    t = dt*[1:size(conWV,1)]*1000;            cmap = jetm(length(PICK_KLUSTERS)); % because
    tbins = [0.5:2:51]; NCOL = 7;
    
    for j=1:length(PICK_KLUSTERS); % plot picked clusters
        i = PICK_KLUSTERS(j);
        ind = find(assign == i);
        waves = conWV(:,ind(1:min(MAXSPIKESPLOT,length(ind))));
        
        % OVERLAYED PLOT
        figure(ifid);
        lh = mplot(t, waves', 'Color', cmap(j,:));
        set(lh, 'ButtonDownFcn', {@raise_me, lh});hold on;
        sleg{j} = ['K: ' num2str(i) ' Ns: ' num2str(K_N(j))];
        
        FID = ifid+1;
        figure(FID); %SUBPLOTS
        h(j) = subplot(length(PICK_KLUSTERS),NCOL,[NCOL*j-(NCOL-1) NCOL*j-5 ]);
        %                 mplot(t, waves', 'Color', cmap(j,:)); hold on;
        [lh,ph] = errorarea(t,mean(waves',1), std(waves',1,1));
        set(lh, 'Color', brighten(cmap(j,:), -0.6), 'ZData', repmat(j, size(get(lh,'XData'))));
        set(ph, 'FaceColor', cmap(j,:), 'ZData', repmat(j, size(get(ph,'XData'))), 'FaceAlpha', 0.8);
        YRANGE(2) = max(YRANGE(2),max(max(waves)));                YRANGE(1) = min(YRANGE(1),min(min(waves))); % set y range to be the same
        
        %%% PLOT AUTOCORR
        subplot(length(PICK_KLUSTERS),NCOL,NCOL*j-4);
        xcorrmk(PICK_KLUSTERS(j),PICK_KLUSTERS(j),stimes,assign,1/dt,'range',40,'fid',FID,'color',cmap(j,:)) % 1ms bin
        
        %%% PLOT ISI
        i = PICK_KLUSTERS(j);
        ind = find(assign == i);
        orderstimes = sort(stimes(ind));
        temp = diff(orderstimes)*dt*1000;
        [a x] = hist(temp,tbins); % ms
        subplot(length(PICK_KLUSTERS),NCOL,NCOL*j-3);
        line(x(1:end-1),a(1:end-1),'color',cmap(j,:)); hold on;
        axis tight; xlabel('ISI (ms)')
        
        %%% PLOT BURST
                subplot(length(PICK_KLUSTERS),NCOL,NCOL*j-2);
                ISIN = 2; % set to 2 for comparing spike isi before and after
               temp2 = diff(orderstimes(ISIN:end))*dt*1000;
        loglog(temp(1:end-(ISIN-1)),temp2,'o'); % plot ISI before vs ISI after each spike
        s = sprintf('Bursts - ISI before vs to %dst spk after ',ISIN-1);
        title(s); xlabel('bef ISI (ms)'); ylabel('aft ISI (ms)')
%               nbins = [10 10]
%         hist3c(temp(1:end-1),temp2, nbins,2)
        xlim([1 10^4])
        ylim([1 10^4])

        % events that occur within ABSREF
        ABSREF = 3;
        line([1 1]*ABSREF,[0 max(a)],'color',[0 0 0]); % indicate absolute refactor values with line
          %%% PLOT SPIKES VIOLATING ABSREF
        subplot(length(PICK_KLUSTERS),NCOL,[NCOL*j-1 NCOL*j]);
        ind = find(assign == i);
        ind0 = ind(find(temp<ABSREF)); NBAD = length(ind0);
        ind1 = ind(find(temp<ABSREF)+1);
        waves = conWV(:,ind);
        % plt events that violate ABSREF
        lh= mplot(t, conWV(:,ind)', 'Color', [0 0 0]); set(lh, 'ButtonDownFcn', {@raise_me, lh});hold on; axis tight % all spikes in cluster
        lh = mplot(t, conWV(:,ind0)', 'Color', [0 1 0]);set(lh, 'ButtonDownFcn', {@raise_me, lh});hold on; axis tight% spike before
        lh = mplot(t, conWV(:,ind1)', 'Color', [1 0 0]');set(lh, 'ButtonDownFcn', {@raise_me, lh});hold on; axis tight% spike after
        title(['NBAD: ' num2str(NBAD) ' frac: ' num2str(NBAD/length(ind)*100) '%'])
        %
        %                 % UNFINISHED (memory runs out) plot klusters accross time
        %                 figure(ifid+1);
        %                 %                              h2(i) = subplot(length(NK),1,i);
        %                 temp(:)=0;
        %                 temp( stimes(ind)) = 1;
        %                 temp = smooth(temp,10e-3/dt);
        %                 plot(stimes(ind)./dt,i*ones(size(ind)),colordot(i,:)); hold on;
        %
    end
    figure(ifid);                       axis tight;    plotset(1);  legend(sleg);plotAnn(sAnn,ifid,2);
    figure(ifid+1);    linkprop(h,{'box','TickDir','Color','XLim','YLim'});      plotset(1); set(gcf,'CurrentAxes',h(1)); ylim(YRANGE);     xlim([t(1) t(end)]);      plotAnn(sAnn,ifid+1,2);
    
end


if bCrossCorr %             cross correlate klusters
    xcorrmk(PICK_KLUSTERS',PICK_KLUSTERS',stimes,assign,1/dt,'range',200,'bsz',5,'fid',100) ; % 5ms bin
    xcorrmk(PICK_KLUSTERS',PICK_KLUSTERS',stimes,assign,1/dt,'range',20,'fid',10) % 1ms bin
end


% CHECK violation of Absolute Refactory period given Poisson process of same (fixed rate)
%             %%% UNFINISHED (ADD to xcorr or figure out how to compare it to integral of
%             %%% autocorrelation
%             ABSREF = 2; % ms;
%             TotalTime = 0; % msec
%             for i = 1:length(datafile) % simulate the number of spikes that should be within ABSREF ms window
%                 TotalTime = TotalTime + datafile(i).MAXTRIGGER *datafile(i).maxsample*datafile(i).dt *1000;% note: total time must be in the same bin size as bins in autocorrelation
%             end
%
%             for i =1:length(NK)
%                 POISS_in_ABSREF(i) = PoissRef(K_rate(i),TotalTime,ABSREF);
%             end