function [usf reSorted spikes] =  sortExtSpikes(extractedSpkFilenames,breSort)
% function [usf reSorted spikes] =  sortExtSpikes(extractedSpkFilenames,SAVEPATH,breSort)
% INPUT:
%      extractedSpkFilenames - cell of file names with extracted data
%      (generated with extractDAQspikes
%      SAVEPATH - optional location for saving extracted spikes
%      breSort - set to 1 to force resort of data even if same data has
%      been sorted before. If set to 0 will reload and plot last sort
%      (assuming that it is at SAVEPATH)
% OUTPUT:
%      usf - savefilename for unit data
%      spikes - struct of chronux klustering (is empty for klustakwik
%      cluster)
%      reSorted - 1 of reSorted
%% Load in spiketimes file
% (generated with extractDAQSpikes.m)
% uses Chrounx soring
% BA
    SORTEDSPIKES_SAVEPATH = []; 

% options
bUSERMANUAL = 1; % User manually defines klusters
bCrossCorr = 0; % cross correlation spike times of units

DFIGPOS = [ 905    43   766   943];
DFIGOFF = [.9 1 1 1];

%% get filename and path info
ANALYZED_DAQSAVEPATH = [];
SORTEDSPIKES_SAVEPATH = []; 
s_RIGSPECIFIC_SPIKESORT_CONFIG; % load defaults

FILENAME = extractedSpkFilenames{1};
SAVEPATH =  SORTEDSPIKES_SAVEPATH;

LOADPATH = fullfile(ANALYZED_DAQSAVEPATH, 'ExtractedSpikes\'); % default save path


if ~exist('breSort','var'); breSort = 0; end; % by default reload last sorted data


SAVESUFFIX = '_SSunitdata';
sf = sprintf('%s_SSFig_%s',FILENAME);
usf = sprintf('%s%s',FILENAME,SAVESUFFIX); % spikesorted save filename


 

%  check files are compatible
for i = 1:length(extractedSpkFilenames)
    % NOTE could save memory by not loading all extSpikes at once but
    % sequentially (then wouldn't predefine conWV size
    lf = fullfile(LOADPATH,extractedSpkFilenames{i});
    if i==1
        load(lf,'extSpikes');
        load(lf,'rawdatainfo');
        
    else % this is silly but don't know how else to make an array of structs that are loaded (without make them substructs of another struc)
        temp = load(lf,'extSpikes');
        extSpikes(i) = temp.extSpikes;
        temp = load(lf,'rawdatainfo');
        rawdatainfo(i) = temp.rawdatainfo;
    end
    
    load(lf,'extSpikesparams');
    if i==1
        dt = rawdatainfo.dt;
        WOI = extSpikesparams.WOI;
        TET = extSpikesparams.TET;
        maxlevel = extSpikesparams.maxlevel;
        bZSCORE = extSpikesparams.bZSCORE;
    else
        
        if dt ~= rawdatainfo(i).dt |...
                WOI ~= extSpikesparams.WOI |...
                TET ~= extSpikesparams.TET |...
                maxlevel ~= extSpikesparams.maxlevel |...
                bZSCORE ~= extSpikesparams.bZSCORE
            ERROR('extractedSpkFilenames do not have the same exSpike Parameters.');
        end
        
    end
    temp =0;
    for j = 1:length(extSpikes(i).spiketimes) % for each electrode
        temp = temp +   size(extSpikes(i).spiketimes(j).ind_sw_trough,1);
    end
    NspikesPerFile(i) = temp;
    
end
clear extSpikesparams;

%% concatinate waveforms,spiketimes
[stimes conWV] = helpconcSpikeData(extSpikes,rawdatainfo); % (will need this even if data already sorted)

sAnn = FILENAME;
bSort = 1;
if ~breSort % load previously sorted data
    if ~isempty(dir(fullfile(SAVEPATH,[usf '.mat'])))
        load(fullfile(SAVEPATH,[usf '.mat']),'spikes','assign');
        bSort = 0;
    end
end
if bSort
    %% dEAL WITH
    % thresholds% *** ASSUME that threshold doesn't change across files
    th = [];
    for i = 1: length(extSpikes)
        for j = 1:length(extSpikes(i).spiketimes)
            th = [th extSpikes(i).spiketimes(j).th];
        end
    end
    
    % keep REAL spiketimes to save later
    for i = 1: length(extSpikes)
        spikesortdata.file(i).spiketimes= nan( NspikesPerFile(i),1);
        spikesortdata.file(i).sweep= nan( NspikesPerFile(i),1);
        spikesortdata.file(i).assign= nan( NspikesPerFile(i),1);
        ind = 1;
        for j = 1:length(extSpikes(i).spiketimes)
            temp = length(extSpikes(i).spiketimes(j).ind_sw_trough);
            spikesortdata.file(i).spiketimes(ind:ind+temp-1)= extSpikes(i).spiketimes(j).ind_sw_trough;
            spikesortdata.file(i).sweep(ind:ind+temp-1)= extSpikes(i).spiketimes(j).sw;
            ind = ind+temp;
        end
    end
    clear extSpikes;
    
    %% deal with outliers using CHX
    clear spikes
    % spikes.waveforms = double((conWV-meanWV)');
    spikes.waveforms = double((conWV)');
    
    spikes.spiketimes = stimes*dt;
    spikes.Fs = 1/dt;
    spikes.threshT= 50; %round(WOI(1)/1e3/dt)+1+50;
    spikes.threshV = [mean(th) inf];
    spikes.thresh = spikes.threshV(1); % this seems to have been added in Dan Hill's version
    
    sortingAlg = 'CHX';
    %********** Chronux sort % see SpikeSortingDemo.m for details
    clear  unit;
    
%     spikes = ss_outliers(spikes);
%     figure(3);
%     subplot(2,1,1);plot(spikes.waveforms'); axis tight; title('Centered Data w/ Outliers Removed');
%     subplot(2,1,2);plot(spikes.outliers.waveforms'); title('Outliers');axis tight;

%     opts.divisions = 32;
%     spikes = ss_kmeans(spikes, opts);
%     opts.divisions = 32;
    spikes = ss_kmeans(spikes);
    figure(6);clf;set(gcf,'Position',DFIGPOS);
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
end

figure(8);clf; set(gcf,'Position',DFIGPOS.*DFIGOFF); colormap jet; % THis plot overlaps with plotKluster but has spiking over time so I kept it
showclust(spikes, assign);

[ifid PICK_KLUSTERS]=local_plotclusters(assign, conWV, stimes,dt,sAnn,bCrossCorr) ;
clear conWV meanWV

if bSort % if sorted then must also save
    
    if ~exist('unit','var'); unitL = 0; else L = length(); end;
    Kluster2extract = PICK_KLUSTERS; % can extract a subset of all clusters
    for i=1:length(Kluster2extract) % extract
        I = find(assign==Kluster2extract(i));
        unit(unitL+i).spiketimes = spikes.spiketimes(I);
        unit(unitL+i).index_in_spikes = I;
    end
    
    spikesortdata.extractedSpkFilenames = extractedSpkFilenames;
    %     spikesortdata.unit = unit; % not sure I need this
        clear  unit;

    OFF = 1;
    for i = 1: length(NspikesPerFile)
        spikesortdata.file(i).assign= assign(OFF:NspikesPerFile(i)+OFF-1);
        OFF = OFF + NspikesPerFile(i);
        % check is same length as spiketimes
        if ~ isequal(size(spikesortdata.file(i).assign),size(spikesortdata.file(i).spiketimes))
            error('spikestimes and assigns are different size, there is an error') % 
            % this shouldn't happen unless spikes are somehow misallocated
            % to the wrong file
        end
    end
    try
    set(0,'CurrentFigure',ifid+1);    savefigure(SAVEPATH,'',sf,'eps');
    catch ME
       getReport(ME);
    end
        save(fullfile(SAVEPATH,usf),'spikesortdata','spikes','assign','rawdatainfo');

end


reSorted = bSort;

    function [ifid PICK_KLUSTERS]=local_plotclusters(assign, conWV, stimes,dt,sAnn,bCrossCorr)
        % plots klusters with greater than X spikes
        % ifid (figure number of plots)
        % PICK_KLUSTERS klusters that are ploted
        DFIGPOS = [ 635    41   686   936];
        DFIGOFF = [.9 1 1 1];
        
        ifid = 4;
        NK = unique(assign); % plot
        NK = NK((NK>0)); % 0 is outliers
        
        % don't bother with klusters with less than 50 spikes
        clear K_N;
        for i= 1:length(NK); K_N(i) = sum(assign == NK(i)); end
        PICK_KLUSTERS = NK((K_N>20));  % only take cluster with more than X spikes
        display('Only take Klusters with > 20 spikes');
        
        if 0 % *** analyze klusters
            figure(ifid); clf;   set(gcf, 'Renderer','painters','Position',DFIGPOS);
            figure(ifid+1); clf;set(gcf, 'Renderer','painters','Position',DFIGPOS.*DFIGOFF.^2);
            %             temp = zeros(1,max(stimes),'int16'); % for creating spikerate as a function of time
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
                ind0 = ind((temp<ABSREF)); NBAD = length(ind0);
                ind1 = ind(find(temp<ABSREF)+1);
                %                 waves = conWV(:,ind);
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
    end
    function [stimes conWV meanWV] =  helpconcSpikeData(extSpikes,rawdatainfo)
        % datafile(N) - struct containing spiketimes and Waveform (WV) (optinal)
        % for N data file in an experiment Data containing
        %
        
        %concatenate time sesimultaneous
        bSUBTRACT_MEAN = 1; % remove common mode noise by subtracting mean of all tetrodewaveforms.
        conWV = [];meanWV = [];stimes = [];
        % TO DO declare all these variables in there full size first.
        
        for m = 1:length(extSpikes)
            for j = 1:length(extSpikes(m).spiketimes)
                
                if nargout >1
                    conWV = [conWV reshape(extSpikes(m).spiketimes(j).WV,size(extSpikes(m).spiketimes(j).WV,1)*size(extSpikes(m).spiketimes(j).WV,2),size(extSpikes(m).spiketimes(j).WV,3))];
                    if bSUBTRACT_MEAN % subtract mean waveform
                        meanWV = [meanWV repmat(squeeze(mean(extSpikes(m).spiketimes(j).WV,2)),size(extSpikes(m).spiketimes(j).WV,2),1)];
                    else meanWV = 0; end
                end
                swlength = rawdatainfo(m).maxsample;
                
                fileOffset = 0;
                for i = m-1:-1:1;            fileOffset = fileOffset+ rawdatainfo(i).maxtrigger*rawdatainfo(i).maxsample;        end % make different files continous by adding time offset before spiketime of file
                stimes = [stimes; extSpikes(m).spiketimes(j).ind_sw_trough+(extSpikes(m).spiketimes(j).sw-1)*swlength+fileOffset]; %
            end
        end
    end
end
