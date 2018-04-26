function [savefilenames breExtracted] = extractDAQSpikes(STF,TET,Extrparams,bNoLoad,SAVEPATH)
% function [savefilenames breExtracted] = extractDAQSpikes(STF,TET,Extrparams,bNoLoad,SAVEPATH))
% reads in 1 or many data files (data in DAQ toolbox .daq format)
% size(TET) channels are read in from each files
%  - data is filtered
% -spiketimes and waveforms extracted
% INPUT:
% STF - struct with fields
%              filename - name of data file to read and use
%              MAXTRIGGER - (optional) scalar with number of tiggers (Starting at 1 to read)
%             size of STF, N  is the number of data files to read in
% TET - vector of chns to read in from datafile
%       waveform of each chn will be extracted if ANY of the chns has a threshold crossing
%       (usually Tetrode is entered but could be any number of channels)
% Extrparams
%     .TH_STD (opt) - threshold is sd for spike extraction (default 4)
%     .invt (opt) - set to -1 to invert trace and looke for upward going spikes
%     .bCLEAN60HZ (opt) default 0
%     .maxlevel (opt) for high pass filter before spike extraction (default 5)
%     .WOI (opt)  %[0.2 0.4]; % size of waveform to extract in ms
% bNoLoad (opt) - Disable loading of previous filtered or extracted data
% SAVEPATH - directory STF files are located
%
% OUTPUT:
%  savefilenames - file name of saved data
%  breExtracted = 1 if data was reExtracted
%
% ALSO  SAVES this struct
% extSpikes - struct
%
%               length(extSpikes) =  1  (only one rawextSpikes perextSpikes)
%               size(extSpikes.spiketimes) = size(TET) %i.e. number ofchannels
% MOST IMPORTANT fields (for each spike)
%         datafile.spiketimes(j).sw  % sweep or trigger of spike
%         datafile.spiketimes(j).ind_sw  % index in spweek if TROUGH in  spike
%         datafile.spiketimes(j).ind  % index in file if sweeps are concatenated
%         datafile.spiketimes(j).ind_sw_trough
%         datafile.spiketimes(j).RANGE % measure of amplitude of spike
%         datafile.spiketimes(j).WV  % waveform
%         datafile.spiketimes(j).th = th;
% (fields that applay to all spikes)
%     datafile.rawdata.dt = dt;
%     datafile.rawdata.filename = filename;
%     datafile.rawdata.swlength = maxsample
%     datafile.rawdata.maxtrigger = maxtrigger;
%     datafile.rawdata.LOADPATH = PATH;
% parameters fields
%     datafile.params.TET = TET;
%     datafile.params.bCLEAN60HZ = bCLEAN60HZ;
%     datafile.params.bZSCORE = bZSCORE;
%     datafile.params.bPLOT = bPLOT;
%     datafile.params.invt = INVt;
%     datafile.params.WOI = WOI; % size
%     datafile.params.TH_STD = TH_STD;
%     datafile.params.WVLT_HPFILTER = HPFILTER;
%     datafile.params.invt = invt ;
%

%default general params
if ~exist('bNoLoad','var'); bNoLoad = 0; end
bPLOT = 0;
bZSCORE = 0;
flagClean60Hz = []; % Set to empty if you want a subset of the data to be tested for 60HZ
bChopUpData = 1;% because arrays get too big and wavefilter chokes
ARBLIMIT = 1*32e3;

DEFLOADPATH = 'C:\Documents and Settings\Bass\My Documents\Scanziani Lab\VC proj\Data\Rawdata';

% set default Extraction parameters
DEFAULTTH_STD = 4;
DEFAULTINVT = 1; % no invert
DEFAULT60HZ = 0;
DEAFAULTWOI = [0.2 0.4]; % size of waveform to extract in ms
DEFAULTmaxlev = 5;
if ~exist('Extrparams','var')
    TH_STD = DEFAULTTH_STD;  % TH_STD = 3; % threshold in STDs
    invt = DEFAULTINVT;
    bCLEAN60HZ = DEFAULT60HZ;
    WOI= DEAFAULTWOI;
    maxlevel = DEFAULTmaxlev ; % cutoff frequency for high pass filter of rawdata = (1/dt)/(2^(maxlevel+1))
else
    if isfield(Extrparams,'TH_STD'); TH_STD = Extrparams.TH_STD;
    else TH_STD = DEFAULTTH_STD; end
    if isfield(Extrparams,'invt'); invt = Extrparams.invt;
    else invt = DEFAULTINVT; end
    if isfield(Extrparams,'bCLEAN60HZ'); bCLEAN60HZ = Extrparams.bCLEAN60HZ;
    else bCLEAN60HZ = DEFAULT60HZ; end
    if isfield(Extrparams,'WOI'); WOI = Extrparams.WOI;
    else WOI = DEAFAULTWOI; end
    if isfield(Extrparams,'maxlevel'); maxlevel = Extrparams.maxlevel;
    else maxlevel = DEFAULTmaxlev; end
end

% set save locations
DAQSAVEPATH = [];s_RIGSPECIFIC_SPIKESORT_CONFIG;
LOADPATH = DAQSAVEPATH;

SAVEPATHFILTERED = fullfile(LOADPATH,'Filtered\');
if isempty(dir(SAVEPATHFILTERED)); mkdir(LOADPATH,'Filtered') ; end
if ~exist('SAVEPATH','var')
    ANALYZED_DAQSAVEPATH = []; s_RIGSPECIFIC_SPIKESORT_CONFIG;
    SAVEPATH = fullfile(ANALYZED_DAQSAVEPATH, 'ExtractedSpikes\'); % default save path
end

% start extracting
for m = 1:length(STF)
    temp = strfind(STF(m).filename,'.');  if isempty(temp) temp = length(STF(m).filename) ; else temp  = temp(end)-1; end % remove file extension
    SAVEFILENAME = sprintf('%s_SpkTimes_TET%s_%sSD%db60%dbZ%d%dWOI%dinvt',STF(m).filename(1:temp(end)),strrep(num2str(TET),'  ','_'),regexprep(num2str(TH_STD),'\.','x'),bCLEAN60HZ,bZSCORE,WOI(1)*100,WOI(2)*100,invt);
    SAVEFULLFILENAME = fullfile(SAVEPATH,SAVEFILENAME);
    savefilenames{m} = SAVEFILENAME;
    clear extSpikes;
    
    bload = 1;
    % check if spikes have already been extracted from this file
    if ~isempty(dir([SAVEFULLFILENAME '.mat']))&& ~bNoLoad 
        try
            load(SAVEFULLFILENAME,'extSpikesparams');
            printf('Loading Previously Extracted Data:\n %s',...
                SAVEFILENAME)
            
            % requested parameters
            params.TET = TET;
            params.bCLEAN60HZ = bCLEAN60HZ;
            params.bZSCORE = bZSCORE;
            params.WOI = WOI;
            params.TH_STD = TH_STD;
            %         params.WVLT_HPFILTER = HPFILTER; % can't check this because requires dt which requires loading the raw data file
            params.maxlevel = maxlevel;
            params.invt = invt;
            % loaded params
            ldparams.TET = extSpikesparams.TET ;
            ldparams.bCLEAN60HZ =   extSpikesparams.bCLEAN60HZ ;
            ldparams.bZSCORE =  extSpikesparams.bZSCORE;
            ldparams.WOI = extSpikesparams.WOI;
            ldparams.TH_STD = extSpikesparams.TH_STD;
            %         ldparams.WVLT_HPFILTER =   extSpikesparams.WVLT_HPFILTER;
            ldparams.maxlevel =   extSpikesparams.maxlevel;
            ldparams.invt =   extSpikesparams.invt;
            
%             disp('Checking Parameters');
            if isempty(comp_struct(params,ldparams,[],[],[],[],0));
                bload = 0;
            else
                warning('Loaded file params do NOT match params. So the way savefilename is names should be changed to fix this.')
            end
        catch ME
            getReport(ME)
        end
    end
    
    
    if bload % reload rawdata if extSpikes doesn't match requested extraction paramters or doesn't exist
        LOADFILE = fullfile(LOADPATH,STF(m).filename);
        rawdatainfo.filename = STF(m).filename;
        if isfield(STF,'MAXTRIGGER')
            if length(STF(m).MAXTRIGGER)==1 && STF(m).MAXTRIGGER>1;TrigRange = [1 STF(m).MAXTRIGGER];
            else TrigRange = STF(m).MAXTRIGGER; end
            [data dt] = loadDAQData(LOADFILE,TET,TrigRange);
        else
            [data dt] = loadDAQData(LOADFILE,TET);
        end
        maxsample = size(data,1); % sweeplength       
        maxtrigger = size(data,2); %update MAXTRIGGER
        
        sAnn = strfind(LOADFILE,'\');sAnn = LOADFILE(sAnn(end)+1:end); % for annotation of figures
        %%
        
        [numpoints, numsweeps, numwires] = size(data);
        f = fullfile(SAVEPATHFILTERED,[STF(m).filename '_TET ' strrep(num2str(TET),'  ','_') 'filtered']);
        if ~isempty(dir([f '.mat']))&& ~bNoLoad % check if these data have been filtered already
            load(f);
            printf('Loading Previously Filtered Data\n Wavelet HIGH-PASS Filter: %s Hz dfilter: %s',num2str(HPFILTER) , num2str(dfilter));
        else
            
            bChopped = 0; % flag set to 1 if data gets chopped
            if bChopUpData % transform into chopped form (later transform back)
                if numpoints >  ARBLIMIT% arb limit
                    if numsweeps ==1 % only reshaping if there is just 1 sweep
                        bChopped =1;
                        for i =1:numwires
                            realnumpoints = numpoints;
                            numsweeps =ceil(numpoints/ARBLIMIT);
                            
                            if i==1; tempd = zeros(ARBLIMIT,numsweeps,numwires,class(data)); end
                            temp = data(:,:,i);
                            temp = [temp(:); temp(end)*ones([numsweeps*ARBLIMIT-numpoints],1)];
                            tempd(:,:,i) = reshape(temp,ARBLIMIT,numsweeps,1);
                        end
                        data = tempd;
                        clear tempd;
                    end
                end
            end
            
            %%% test if need to remove 60Hz
            if isempty(flagClean60Hz)||flagClean60Hz==1
                prepdata(squeeze(data(1:1/dt,1,1))',1,dt,[NaN,NaN,60 1]);
                flagClean60Hz = 0;
                if 0 % if you want a dialog to popup each time and ask whether to clean 60Hz
                    button = questdlg('Clean 60Hz?','Filtering',...
                        'Yes','No to All','Yes to All','No to All');
                    
                    if strcmp(button,'Yes')
                        flagClean60Hz = 1;
                    elseif strcmp(button,'No to All')
                        flagClean60Hz = 0;
                    elseif strcmp(button,'Yes to All')
                        flagClean60Hz = 2;
                    end
                end
                if flagClean60Hz; bCLEAN60HZ=1; else bCLEAN60HZ=0; end
            end
            
            
            if bCLEAN60HZ
                dfilter = [NaN,NaN,60 1]; clear df;
                tic
                display('Cleaning 60Hz')
                for i =1 : size(data,3)
                    data(:,:,i) = prepdata(squeeze(data(:,:,i))',0,dt,dfilter)';
                end
                toc
            else
                dfilter = [NaN,NaN,60 0];
            end
            
            HPFILTER = (1/dt)/(2^(maxlevel+1));
            display(['Wavelet HIGH-PASS Filter: ' num2str(HPFILTER) ' Hz'])
            tic;df = wavefilter(data(:,:,:),maxlevel);toc
            
            if bChopped % transform back to correct form if chopped up
                for i =1:numwires
                    if i==1; tempd = zeros(realnumpoints,1,numwires,class(data)); end
                    temp = df(:,:,i);
                    tempd(:,1,i) = temp(1:realnumpoints);
                end
                df = tempd;
                clear tempd;
            end
            
            if size(data,2)==1 % for case of only 1 file make 3D to be compatible
                df  = permute(df, [1,3,2] );
            end
            save(f,'df','dfilter','HPFILTER');
        end
        for i = 1:size(df,3)  % DATA?
            STDEV(i) = median(std(squeeze(df(:,:,i))));
            % BA 091909 using median rather than mean is more robust to outliers
            % **however if there is only 1 sweep is still not robust Median
            % absolute deviation would be better, only slight disadvantage
            % of this is that you can't specify threshold in STDs but must use MADs
            % instead
            if bZSCORE
                df(:,:,i) = df(:,:,i)/STDEV(i);
                data(:,:,i) = data(:,:,i)/STDEV(i);
            end
        end% convert to Z-score
        
        if bPLOT
            clear h;
            ifid = 1;
            figure(ifid);clf; j =1;
            for i =1: length(TET); h(i) =subplot(length(TET),1,i); plotdata(data(:,j,i),dt,'color','b');plotdata(df(:,j,i),dt,'color','r');
                ylabel(['SD = ' num2str(STDEV(i)*1000,'%1.0f') 'uV' ]); axis tight
            end; linkaxes(h,'x');    linkprop(h,{'box','TickDir','Color'});
            plotset(1); plotAnn(sAnn,ifid,3);
        end% % plot tetrode
        
        clear data;
        df = df*invt;
        
        for i = 1:size(df,3) % THRESHOLD for each site
            if bZSCORE
                [extSpikes.spiketimes(i).th extSpikes.spiketimes(i).ind_sw extSpikes.spiketimes(i).sw]=threshdata(squeeze(df(:,:,i)), dt,-TH_STD,[],1); %% asking for threshold input
            else
                [extSpikes.spiketimes(i).th extSpikes.spiketimes(i).ind_sw extSpikes.spiketimes(i).sw]=threshdata(squeeze(df(:,:,i)), dt,-TH_STD*STDEV(i),[],1); %% asking for threshold input
            end
        end
        for i = 1:size(df,3) % spikes with a WOI around them
            temp = [extSpikes.spiketimes(i).ind_sw-round(WOI(1)/1e3/dt)+1 extSpikes.spiketimes(i).ind_sw+round(WOI(2)/1e3/dt)+1];
            includeind = find(temp(:,1)>=1&(temp(:,2)<size(df,1)));
            extSpikes.spiketimes(i).ind_sw = extSpikes.spiketimes(i).ind_sw(includeind);
            extSpikes.spiketimes(i).sw = extSpikes.spiketimes(i).sw(includeind);
        end
        for i = 1:size(df,3) %% estimate amplitude of spike (used to determine which spike to keep when there is a collision)
            for j=1:length(extSpikes.spiketimes(i).ind_sw)
                extSpikes.spiketimes(i).RANGE(j) = range(df(extSpikes.spiketimes(i).ind_sw(j):extSpikes.spiketimes(i).ind_sw(j)+round(0.5/1e3/dt)-1,extSpikes.spiketimes(i).sw(j),i));
            end
        end
        
        MINWINDOW = (round(sum(WOI)/1e3/dt) +1); % min time between events
        for k = 1:3 % need to figure out why this needs to be run more than once.
            for i = 1:size(df,3) % find simultaneous events on same electrode (take larger one)
                if ~isempty(extSpikes.spiketimes(i).sw) % so doesn't crash if no events are found
                    extSpikes.spiketimes(i).ind = extSpikes.spiketimes(i).ind_sw+(extSpikes.spiketimes(i).sw-1)*size(df,1);
                    mask = [1; (diff(extSpikes.spiketimes(i).ind)> MINWINDOW)];
                    ind_sim = find(mask==0);
                    if ~isempty(ind_sim)
                        for j =1:length(ind_sim)
                            if extSpikes.spiketimes(i).RANGE(ind_sim(j)-1)<extSpikes.spiketimes(i).RANGE(ind_sim(j))
                                mask(ind_sim(j)-1) = 0; mask(ind_sim(j)) = 1; end % keep largest event
                        end
                    end
                    extSpikes.spiketimes(i).sw = extSpikes.spiketimes(i).sw(mask==1); % REMOVE spikes within WOI of each other (in same electrode)
                    extSpikes.spiketimes(i).ind_sw = extSpikes.spiketimes(i).ind_sw(mask==1) ;
                    extSpikes.spiketimes(i).ind = extSpikes.spiketimes(i).ind(mask==1) ;
                    extSpikes.spiketimes(i).RANGE = extSpikes.spiketimes(i).RANGE(mask==1) ;
                end
            end
        end
        
        % **  NOT TESTED
        % find events detected on multiple sites
        for i = 1:size(df,3); extSpikes.spiketimes(i).simult = []; end
        
        for j = 1:  maxtrigger % for each sweep
            for i = 1:size(df,3) % for each electrode site
                swevent = find (extSpikes.spiketimes(i).sw==j);
                for k = 1: length(swevent); % for each spike
                    eventTime = extSpikes.spiketimes(i).ind_sw(swevent(k)); % time of a spike in a sweep
                    for l = i+1:size(df,3) % for all OTHER sites
                        swevent2 = find(extSpikes.spiketimes(l).sw==j);
                        eventTime2 = extSpikes.spiketimes(l).ind_sw(swevent2); % time of a spikeS in a sweep
                        temp = swevent2((abs(eventTime2-eventTime)<MINWINDOW));
                        if ~isempty(temp)
                            for n =1: length(temp)
                                extSpikes.spiketimes(i).simult = [extSpikes.spiketimes(i).simult; swevent(k) temp(n) l extSpikes.spiketimes(i).RANGE(swevent(k)) extSpikes.spiketimes(l).RANGE(temp(n)) ]; % index of simultaneous event: THIS site, OTHER site, OTHER site #
                            end
                        end
                    end
                end
            end
        end
        % remove events detected on multiple sites
        clear temp;for j = 1:size(df,3);   temp(j).mask = ones(length(extSpikes.spiketimes(j).sw),1);end % initialize masks
        for j = 1:size(df,3)
            if ~isempty(extSpikes.spiketimes(j).simult)
                temp_sw = unique(extSpikes.spiketimes(j).simult(:,1));
                for k = 1 : length(temp_sw)% all simultaneous events
                    temp_i = find(extSpikes.spiketimes(j).simult(:,1)==temp_sw(k));
                    MAX2 = max(extSpikes.spiketimes(j).simult(temp_i,5));
                    if MAX2 < extSpikes.spiketimes(j).simult(temp_i(1),4); % if max of OTHER sites is less than current site remove event from other sites
                        for i =1:length(temp_i); temp(extSpikes.spiketimes(j).simult(temp_i(i),3)).mask(extSpikes.spiketimes(j).simult(temp_i(i),2)) = 0; end
                    else
                        temp(j).mask(temp_sw(k)) = 0; % remove event from current site
                        if temp_i>1
                            temp2 = extSpikes.spiketimes(j).simult(temp_i(1)-1+find(extSpikes.spiketimes(j).simult(temp_i,5)~=MAX2),3); % remove from OTHER sites
                            for i = 1:length(temp2); temp(extSpikes.spiketimes(j).simult(temp2(i),3)).mask(extSpikes.spiketimes(j).simult(temp2(i),2)) = 0;   end
                        end
                    end
                end
            end
            sum(temp(j).mask)
        end % set mask to zero for simultaneous spikes on sites where the amplitude is not the largest
        for j = 1:size(df,3)% REMOVE spikes within WOI of each other (in same electrode)
            extSpikes.spiketimes(j).sw = extSpikes.spiketimes(j).sw(temp(j).mask==1);
            extSpikes.spiketimes(j).ind_sw = extSpikes.spiketimes(j).ind_sw(temp(j).mask==1) ;
            extSpikes.spiketimes(j).ind = extSpikes.spiketimes(j).ind(temp(j).mask==1) ;
            extSpikes.spiketimes(j).RANGE = extSpikes.spiketimes(j).RANGE(temp(j).mask==1) ;
        end
        
        %%trouble shooting.. find out events that are less than minwindow apart
        %%(from same electrode)
        % for i=1:200
        %     e = [];
        %     for j = 1:size(df,3)
        %         tempI = extSpikes.spiketimes(j).sw==i;
        %         e = [e ;extSpikes.spiketimes(j).ind(tempI)];
        %     end
        % min(diff(sort(e))) < MINWINDOW
        % pause;
        % end
        
        % check for simulaneous events
        eventTimes = []; 
        for j =1:size(df,3)
            eventTimes = [eventTimes; extSpikes.spiketimes(j).ind];
        end % concatenate events from different sites
        if  min(diff(sort(eventTimes))) < MINWINDOW % check that no simultaneous events till exist
            error('simultaneous events not correctly removed')
        end
        for j = 1: size(df,3)
            extSpikes.spiketimes(j).ind_sw_trough = nan(length(extSpikes.spiketimes(j).ind_sw),1);
            for i = 1:length(extSpikes.spiketimes(j).ind_sw) % get ind of trough for each spike
                wave = df(extSpikes.spiketimes(j).ind_sw(i):extSpikes.spiketimes(j).ind_sw(i)+round(WOI(2)/1e3/dt), extSpikes.spiketimes(j).sw(i),j);
                extSpikes.spiketimes(j).ind_sw_trough(i) = find(wave==min(wave),1,'first')+extSpikes.spiketimes(j).ind_sw(i)-1;
            end
        end% find time of trough for each threshold crossing
        % using this rather then center of mass or something
        % OK decided on event times now extract waveforms
        
        for i = 1:size(df,3) % spikes with a WOI around ind_sw_trough ** SHOULDNT havta do a similar step above and here again
            temp = [extSpikes.spiketimes(i).ind_sw_trough-round(WOI(1)/1e3/dt)+1 extSpikes.spiketimes(i).ind_sw_trough+round(WOI(2)/1e3/dt)+1];
            includeind = find(temp(:,1)>=2&(temp(:,2)<size(df,1)));
            extSpikes.spiketimes(i).ind_sw = extSpikes.spiketimes(i).ind_sw(includeind);
            extSpikes.spiketimes(i).sw = extSpikes.spiketimes(i).sw(includeind);
            extSpikes.spiketimes(i).RANGE = extSpikes.spiketimes(i).RANGE(includeind);
            extSpikes.spiketimes(i).ind= extSpikes.spiketimes(i).ind(includeind);
            extSpikes.spiketimes(i).ind_sw_trough= extSpikes.spiketimes(i).ind_sw_trough(includeind);
        end
        
        try
            for k = 1:size(df,3) % get waveforms for each spike trough
                extSpikes.spiketimes(k).WV = [];
                for i=1:length(extSpikes.spiketimes(k).ind_sw_trough)
                    for j = 1:size(df,3)
                        extSpikes.spiketimes(k).WV(:,j,i)=df(extSpikes.spiketimes(k).ind_sw_trough(i)-round(WOI(1)/1e3/dt):extSpikes.spiketimes(k).ind_sw_trough(i)+round(WOI(2)/1e3/dt),extSpikes.spiketimes(k).sw(i),j);
                    end
                end
            end
        catch ME
            getReport(ME)
            k
            keyboard;
        end
        % save
        rawdatainfo.dt = dt;
        rawdatainfo.swlength = maxsample;
        rawdatainfo.maxsample = maxsample;
        rawdatainfo.maxtrigger = maxtrigger;
        extSpikesparams.TET = TET;
        extSpikesparams.bCLEAN60HZ = bCLEAN60HZ;
        extSpikesparams.bZSCORE = bZSCORE;
        extSpikesparams.bPLOT = bPLOT;
        extSpikesparams.WOI = WOI;
        extSpikesparams.TH_STD = TH_STD;
        extSpikesparams.STD =  STDEV;
        extSpikesparams.WVLT_HPFILTER = HPFILTER;
        extSpikesparams.maxlevel = maxlevel;
        extSpikesparams.invt = invt ;
        
        clear df;
        
        
        % so that extSpikes doesn't have to be loaded in entirity to check if it is the right file
        %     creat another struct with params
        save(SAVEFULLFILENAME,'rawdatainfo','extSpikes','extSpikesparams')
    end
end


breExtracted = bload;

% catch
%
%     keyboard
%
% end