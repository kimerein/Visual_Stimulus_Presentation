function [SAVEFILENAME datafile] = extractDAQ_Spikes(STF,TET,PATH,TH_STD,bNoLoad,INVt,bCLEAN60HZ)
% function [SAVEFILENAME datafile] = extractDAQ_Spikes(STF,TET,PATH,TH_STD,bNoLoad,INVt,bCLEAN60HZ)
% reads in 1 or many data files (data in DAQ toolbox .daq format)
% size(TET) channels are read in from each files
% data is filtered
% spiketimes and waveforms extracted
% Data output as struct
% optionally saves data (not SAVEFILENAMEHEADER is changed from input value to
% include threshold)
%
%
% Depends on: loadDAQData.m
% Use in s_ABC_analyzeUnits script
% INPUT:
% STF - struct with fields
%              filename - name of data file to read and use
%              MAXTRIGGER - (optional) scalar with number of tiggers (Starting at 1 to read)
%             size of STF, N  is the number of data files to read in
% TET - vector of chns to read in from datafile
%       waveform of each chn will be extracted if ANY of the chns has a threshold crossing
%       (usually Tetrode is entered but could be any number of channels)
% PATH - directory STF files are located
% SAVEFILENAMEHEADER - (optional) if exists data will be saved to the directory\filename specified
% TH_STD (opt) - threshold is sd for spike extraction (default 4)
% bNoLoad (opt) - Disable loading of previous filtered or extracted data
% INVt (opt) - set to -1 to invert trace and looke for upward going spikes
%
% OUTPUT:
% datafile - struct
%
%               size(datafile) = size(STF); % i.e. number of files
%               size(datafile.spiketimes) = size(TET) %i.e. number ofchannels
% MOST IMPORTANT fields (for each spike)
%         datafile(m).spiketimes(j).sw  % sweep or trigger of spike
%         datafile(m).spiketimes(j).ind_sw  % index in spweek if TROUGH in  spike
%         datafile(m).spiketimes(j).ind  % index in file if sweeps are concatenated
%         datafile(m).spiketimes(j).RANGE % measure of amplitude of spike
%         datafile(m).spiketimes(j).WV  % waveform
% (fields that applay to all spikes)
%     datafile(m).spiketimes(1).dt = dt;
%     datafile(m).dt = dt;
%     datafile(m).spiketimes(1).swlength = maxsample; % for backward compatiblity put maxsample in 2 places
%     datafile(m).maxsample = maxsample;
%     datafile(m).MAXTRIGGER = maxtrigger;
%     datafile(m).params.TET = TET;
% parameters fields
%     datafile(m).params.bCLEAN60HZ = bCLEAN60HZ;
%     datafile(m).params.bZSCORE = bZSCORE;
%     datafile(m).params.bPLOT = bPLOT;
%     datafile(m).params.WOI = WOI; % size
%     datafile(m).params.TH_STD = TH_STD;
%     datafile(m).params.WVLT_HPFILTER = HPFILTER;
%     datafile(m).LOADPATH = PATH;

% set default parameters
if ~exist('TH_STD'); TH_STD = 4; end % TH_STD = 3; % threshold in STDs
if ~exist('bNoLoad'); bNoLoad = 0; end
if ~exist('INVt'); INVt = 1; end % invert 
flagClean60Hz = []; % set when user decides if cleaning is needed
if ~exist('bCLEAN60HZ'); bCLEAN60HZ = 0;end
bZSCORE = 0;
bPLOT = 0;
WOI= [0.2 0.4]; % size of waveform to extract in ms

maxlevel = 5 ; % cutoff frequency = (1/dt)/(2^(maxlevel+1))

% make savefilename
i = 1;
while i <=length(STF); % get filename(s) and concatanate for savefile
    temp = findstr(STF(i).filename,'_');
    if i ==1; stemp = [STF(i).filename(1:temp(end))];end
    
    stemp = [stemp 'f' STF(i).filename(temp(end)+1:end)];% filename
    if i == 18  % impose a maximum length on filename
        i = length(STF)-2 ;
    end
    i = i+1;
end
SAVEPATH = fullfile(PATH,'Analyzed');
SAVEFILENAME = sprintf('%s%s_SpkTimes_TET%s_%sSD%db60%dbZ%d%dWOI',SAVEPATH,stemp,strrep(num2str(TET),'  ','_'),num2str(TH_STD),bCLEAN60HZ,bZSCORE,WOI(1)*100,WOI(2)*100);
% for backwards compatiblity must be added only when inverted
if INVt==-1; SAVEFILENAME = sprintf('%sINVt',SAVEFILENAME); end

SAVEPATHFILTERED = fullfile(PATH,'Analyzed','Filtered');
% SAVEFILENAME = sprintf('%s_spiketimes_%dSD',SAVEFILENAMEHEADER,TH_STD);

bChopUpData = 1;% because arrays get too big and wavefilter chokes
ARBLIMIT = 1*32e3;


% try
bload = 1;
if ~isempty(dir([SAVEFILENAME '.mat']))& ~bNoLoad
    try
    load(SAVEFILENAME,'datafile');
    display(['Loading Previously Extracted Data\n'...
        SAVEFILENAME])  
    
    % requested parameters
        params.TET = TET;
        params.bCLEAN60HZ = bCLEAN60HZ;
        params.bZSCORE = bZSCORE;
        params.WOI = WOI;
        params.TH_STD = TH_STD;
        %         params.WVLT_HPFILTER = HPFILTER; % can't check this because requires dt which requires loading the raw data file
        params.maxlevel = maxlevel;
        params.INVt = INVt;
        % loaded params
        m=1;
        ldparams.TET = datafile(m).params.TET ;
        ldparams.bCLEAN60HZ =   datafile(m).params.bCLEAN60HZ ;
        ldparams.bZSCORE =  datafile(m).params.bZSCORE;
        ldparams.WOI = datafile(m).params.WOI;
        ldparams.TH_STD = datafile(m).params.TH_STD;
        %         ldparams.WVLT_HPFILTER =   datafile(m).params.WVLT_HPFILTER;
        ldparams.maxlevel =   datafile(m).params.maxlevel;
        ldparams.INVt =   datafile(m).params.INVt;
        
        temp = 'Checking Parameters';
        display(temp);
        if isempty(comp_struct(params,ldparams))
            bload = 0;
        else
            warning('Loaded file params do NOT match params. So the way savefilename is names should be changed to fix this.')
        end
    catch ME
        getReport(ME)
    end
end

if bload
    for m = 1:length(STF)
        LOADFILE = [PATH STF(m).filename];
        datafile(m).filename = LOADFILE;
        if isfield('STF','MAXTRIGGER')
            if length(STF(m).MAXTRIGGER)==1 && STF(m).MAXTRIGGER>1;TrigRange = [1 STF(m).MAXTRIGGER];
            else TrigRange = STF(m).MAXTRIGGER; end
            [data dt] = loadDAQData(LOADFILE,TET,TrigRange);
        else
            [data dt] = loadDAQData(LOADFILE,TET);
        end
        data = data*INVt;
        maxsample = size(data,1); % sweeplength
        
        maxtrigger = size(data,2); %update MAXTRIGGER
        
        sAnn = strfind(LOADFILE,'\');sAnn = LOADFILE(sAnn(end)+1:end);
        %%
        
        [numpoints, numsweeps, numwires] = size(data);
        f = [SAVEPATHFILTERED STF(m).filename '_TET ' strrep(num2str(TET),'  ','_') 'filtered'];
        if ~isempty(dir([f '.mat']))& ~bNoLoad % check if these data have been filtered already
            load(f);
            display(['Loading Previously Filtered Data\n'...
                'Wavelet HIGH-PASS Filter: ' num2str(HPFILTER) ' Hz\n'...
                'dfilter: ' num2str(dfilter)])
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
                
                button = questdlg('Clean 60Hz?','Filtering',...
                    'Yes','No to All','Yes to All','No to All');
                
                if strcmp(button,'Yes')
                    flagClean60Hz = 1;
                elseif strcmp(button,'No to All')
                    flagClean60Hz = 0;
                elseif strcmp(button,'Yes to All')
                    flagClean60Hz = 2;
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
        
        for i = 1:size(df,3) % THRESHOLD for each site
            if bZSCORE
                %             [datafile(m).datafile(m).spiketimes(i).th datafile(m).spiketimes(i).ind_sw datafile(m).spiketimes(i).sw]=threshdata(squeeze(df(:,:,i)), dt,-TH_STD); %% asking for threshold input
                [datafile(m).spiketimes(i).th datafile(m).spiketimes(i).ind_sw datafile(m).spiketimes(i).sw]=threshdata(squeeze(df(:,:,i)), dt,-TH_STD,[],1); %% asking for threshold input
            else
                [datafile(m).spiketimes(i).th datafile(m).spiketimes(i).ind_sw datafile(m).spiketimes(i).sw]=threshdata(squeeze(df(:,:,i)), dt,-TH_STD*STDEV(i),[],1); %% asking for threshold input
            end
        end
        for i = 1:size(df,3) % spikes with a WOI around them
            temp = [datafile(m).spiketimes(i).ind_sw-round(WOI(1)/1e3/dt)+1 datafile(m).spiketimes(i).ind_sw+round(WOI(2)/1e3/dt)+1];
            includeind = find(temp(:,1)>=1&(temp(:,2)<size(df,1)));
            datafile(m).spiketimes(i).ind_sw = datafile(m).spiketimes(i).ind_sw(includeind);
            datafile(m).spiketimes(i).sw = datafile(m).spiketimes(i).sw(includeind);
        end
        for i = 1:size(df,3) %% estimate amplitude of spike (used to determine which spike to keep when there is a collision)
            for j=1:length(datafile(m).spiketimes(i).ind_sw)
                datafile(m).spiketimes(i).RANGE(j) = range(df(datafile(m).spiketimes(i).ind_sw(j):datafile(m).spiketimes(i).ind_sw(j)+round(0.5/1e3/dt)-1,datafile(m).spiketimes(i).sw(j),i));
            end
        end
    
        MINWINDOW = (round(sum(WOI)/1e3/dt) +1); % min time between events
        for j = 1:3 % need to figure out why this needs to be run more than once.
            for i = 1:size(df,3) % find simultaneous events on same electrode (take larger one)
                if ~isempty(datafile(m).spiketimes(i).sw) % so doesn't crash if no events are found
                    datafile(m).spiketimes(i).ind = datafile(m).spiketimes(i).ind_sw+(datafile(m).spiketimes(i).sw-1)*size(df,1);
                    mask = [1; (diff(datafile(m).spiketimes(i).ind)> MINWINDOW)];
                    ind_sim = find(mask==0);
                    if ~isempty(ind_sim)
                        for j =1:length(ind_sim)
                            if datafile(m).spiketimes(i).RANGE(ind_sim(j)-1)<datafile(m).spiketimes(i).RANGE(ind_sim(j))
                                mask(ind_sim(j)-1) = 0; mask(ind_sim(j)) = 1; end % keep largest event
                        end
                    end
                    datafile(m).spiketimes(i).sw = datafile(m).spiketimes(i).sw(find(mask)); % REMOVE spikes within WOI of each other (in same electrode)
                    datafile(m).spiketimes(i).ind_sw = datafile(m).spiketimes(i).ind_sw(find(mask)) ;
                    datafile(m).spiketimes(i).ind = datafile(m).spiketimes(i).ind(find(mask)) ;
                    datafile(m).spiketimes(i).RANGE = datafile(m).spiketimes(i).RANGE(find(mask)) ;
                end
            end
        end
        
        % **  NOT TESTED
        % find events detected on multiple sites
        for i = 1:size(df,3); datafile(m).spiketimes(i).simult = []; end
        
        for j = 1:  maxtrigger % for each sweep
            for i = 1:size(df,3) % for each electrode site
                swevent = find (datafile(m).spiketimes(i).sw==j);
                for k = 1: length(swevent); % for each spike
                    eventTime = datafile(m).spiketimes(i).ind_sw(swevent(k)); % time of a spike in a sweep
                    for l = i+1:size(df,3) % for all OTHER sites
                        swevent2 = find(datafile(m).spiketimes(l).sw==j);
                        eventTime2 = datafile(m).spiketimes(l).ind_sw(swevent2); % time of a spikeS in a sweep
                        temp = swevent2(find(abs(eventTime2-eventTime)<MINWINDOW));
                        if ~isempty(temp)
                            for n =1: length(temp)
                                datafile(m).spiketimes(i).simult = [datafile(m).spiketimes(i).simult; swevent(k) temp(n) l datafile(m).spiketimes(i).RANGE(swevent(k)) datafile(m).spiketimes(l).RANGE(temp(n)) ]; % index of simultaneous event: THIS site, OTHER site, OTHER site #
                            end
                        end
                    end
                end
            end
        end
        % remove events detected on multiple sites
        clear temp;for j = 1:size(df,3);   temp(j).mask = ones(length(datafile(m).spiketimes(j).sw),1);end % initialize masks
        for j = 1:size(df,3)
            if ~isempty(datafile(m).spiketimes(j).simult)
                temp_sw = unique(datafile(m).spiketimes(j).simult(:,1));
                for k = 1 : length(temp_sw)% all simultaneous events
                    temp_i = find(datafile(m).spiketimes(j).simult(:,1)==temp_sw(k));
                    MAX2 = max(datafile(m).spiketimes(j).simult(temp_i,5));
                    if MAX2 < datafile(m).spiketimes(j).simult(temp_i(1),4); % if max of OTHER sites is less than current site remove event from other sites
                        for i =1:length(temp_i); temp(datafile(m).spiketimes(j).simult(temp_i(i),3)).mask(datafile(m).spiketimes(j).simult(temp_i(i),2)) = 0; end
                    else
                        temp(j).mask(temp_sw(k)) = 0; % remove event from current site
                        if temp_i>1
                            temp2 = datafile(m).spiketimes(j).simult(temp_i(1)-1+find(datafile(m).spiketimes(j).simult(temp_i,5)~=MAX2),3); % remove from OTHER sites
                            for i = 1:length(temp2); temp(datafile(m).spiketimes(j).simult(temp2(i),3)).mask(datafile(m).spiketimes(j).simult(temp2(i),2)) = 0;   end
                        end
                    end
                end
            end
            sum(temp(j).mask)
        end % set mask to zero for simultaneous spikes on sites where the amplitude is not the largest
        for j = 1:size(df,3)% REMOVE spikes within WOI of each other (in same electrode)
            datafile(m).spiketimes(j).sw = datafile(m).spiketimes(j).sw(find(temp(j).mask));
            datafile(m).spiketimes(j).ind_sw = datafile(m).spiketimes(j).ind_sw(find(temp(j).mask)) ;
            datafile(m).spiketimes(j).ind = datafile(m).spiketimes(j).ind(find(temp(j).mask)) ;
            datafile(m).spiketimes(j).RANGE = datafile(m).spiketimes(j).RANGE(find(temp(j).mask)) ;
        end
        
        %%trouble shooting.. find out events that are less than minwindow apart
        %%(from same electrode)
        % for i=1:200
        %     e = [];
        %     for j = 1:size(df,3)
        %         tempI = datafile(m).spiketimes(j).sw==i;
        %         e = [e ;datafile(m).spiketimes(j).ind(tempI)];
        %     end
        % min(diff(sort(e))) < MINWINDOW
        % pause;
        % end
        
        % check for simulaneous events
        eventTimes = []; for j =1:size(df,3)
            eventTimes = [eventTimes; datafile(m).spiketimes(j).ind];
        end % concatenate events from different sites
        if  min(diff(sort(eventTimes))) < MINWINDOW % check that no simultaneous events till exist
            error('simultaneous events not correctly removed')
        end
        for j = 1: size(df,3)
            datafile(m).spiketimes(j).ind_sw_trough = nan(length(datafile(m).spiketimes(j).ind_sw),1);
            for i = 1:length(datafile(m).spiketimes(j).ind_sw) % get ind of trough for each spike
                wave = df(datafile(m).spiketimes(j).ind_sw(i):datafile(m).spiketimes(j).ind_sw(i)+round(WOI(2)/1e3/dt), datafile(m).spiketimes(j).sw(i),j);
                datafile(m).spiketimes(j).ind_sw_trough(i) = find(wave==min(wave),1,'first')+datafile(m).spiketimes(j).ind_sw(i)-1;
            end
        end% find time of trough for each threshold crossing
        % using this rather then center of mass or something
        % OK decided on event times now extract waveforms
        
        for i = 1:size(df,3) % spikes with a WOI around ind_sw_trough ** SHOULDNT havta do a similar step above and here again
            temp = [datafile(m).spiketimes(i).ind_sw_trough-round(WOI(1)/1e3/dt)+1 datafile(m).spiketimes(i).ind_sw_trough+round(WOI(2)/1e3/dt)+1];
            includeind = find(temp(:,1)>=1&(temp(:,2)<size(df,1)));
            datafile(m).spiketimes(i).ind_sw = datafile(m).spiketimes(i).ind_sw(includeind);
            datafile(m).spiketimes(i).sw = datafile(m).spiketimes(i).sw(includeind);
            datafile(m).spiketimes(i).RANGE = datafile(m).spiketimes(i).RANGE(includeind);
            datafile(m).spiketimes(i).ind= datafile(m).spiketimes(i).ind(includeind);
            datafile(m).spiketimes(i).ind_sw_trough= datafile(m).spiketimes(i).ind_sw_trough(includeind);
        end

        try
        for k = 1:size(df,3) % get waveforms for each spike trough
            datafile(m).spiketimes(k).WV = [];
            for i=1:length(datafile(m).spiketimes(k).ind_sw_trough)
                for j = 1:size(df,3)
                    datafile(m).spiketimes(k).WV(:,j,i)=df(datafile(m).spiketimes(k).ind_sw_trough(i)-round(WOI(1)/1e3/dt):datafile(m).spiketimes(k).ind_sw_trough(i)+round(WOI(2)/1e3/dt),datafile(m).spiketimes(k).sw(i),j);
                end
            end
        end
        catch
            k
        end
        % save
        datafile(m).spiketimes(1).dt = dt;
        datafile(m).dt = dt;
        datafile(m).spiketimes(1).swlength = maxsample;
        datafile(m).maxsample = maxsample;
        datafile(m).MAXTRIGGER = maxtrigger;
        datafile(m).params.TET = TET;
        datafile(m).params.bCLEAN60HZ = bCLEAN60HZ;
        datafile(m).params.bZSCORE = bZSCORE;
        datafile(m).params.bPLOT = bPLOT;
        datafile(m).params.WOI = WOI;
        datafile(m).params.TH_STD = TH_STD;
        datafile(m).params.STD =  STDEV;
        datafile(m).params.WVLT_HPFILTER = HPFILTER;
        datafile(m).params.maxlevel = maxlevel;
        datafile(m).params.INVt = INVt ;
        
        datafile(m).LOADPATH = PATH;
        clear df;
    end
    
    save(SAVEFILENAME,'datafile')
    
end


% catch
%
%     keyboard
%
% end