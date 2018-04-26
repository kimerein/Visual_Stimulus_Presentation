function analyzeUnitTuningNEW(spkSortFN,STF,UNIT,AnalType,cond,includesweeps,analyzewindow,label,bsave)
% function analyzeUnitTuningNEW(spkSortFN,STF,DAQchns,exptparam,UNIT,AnalType)
% saves unit info and figures
%  INPUT
%     includesweeps -  cell vector length N, where N is the number of rawdaq files in STF
%                      includesweeps{i} = [] % ALL sweeps are included from
%                      file
%                      includesweeps{i} = 0 % No sweeps are included from file
%     analyzewindow - interval in sweep to be analyzed as the stimulus period (in sec)
%               e.g. [0 2]
%          when analyzewindow = []; then the stimulus duration from the PsychStimController will be used.
%     label - optional string to describe analysis it will be saved to unitsummaryCond
%     bsave - saves 2 structs with unit analysis data (see code for more details) unitsummaryCond and unitsummary
% ** BUG with this replaces entire unit table some times don't know why yet
%   AnalType - what type of analysis.
%  'Orientation' will load the files in STF where only Orientation is varied.
%     if "cond" is specified all files like the first file (in STF) that
%     meets cond and changes in Orientation. will be analyzed.
%    to match.
%    only files that match in terms of the StimulusName, Duration and
%    baseline duration (in addition to any cond specified) are considered
%   cond(n) is a struct that defines what conditions should be analyzed
%     it contains fields for different condtions:
%     {'.bMask','.Mask','.contrast','.spfreq','.tempFreq','.squaregratings'};conditions can be added in future
%     0) if a condition is specified only that condition will be used
%     1) any conditions not specified, will be defined by the first STF file that meets the AnalType critera (Templatefile)
%     2) cond fields can be to 'all' this will include all conditions of
%     this kind
%                e.g. cond(1).mask=1; cond(2).contrast = 1; cond(1).sqgrating = 0;
%                cond(2).mask=0; cond(3).contrast = 1; cond(1).sqgrating = 0;
%                cond.mask = 'Gau'; %'Inv' %'Ape' % 'Non'
% ********** ONLY TESTED with mask so far
%
%%% BUG when STF is not same length as rawdatainfo causes problem
% TO DO:a
% saving unit needs work (and is not updated if spikesorts are removed)


if ~exist('bsave','var'); bsave = 1; printf('*****************************\n \t\t\tUnit data will be SAVED');  end
if ~exist('includesweeps','var') || isempty(includesweeps); includesweeps = cell(size(STF)); end
if ~exist('analyzewindow','var'); analyzewindow = []; end
if ~exist('label','var'); label = []; end


% TO DO: if spkSortFNnot specified use newest sort
TABLES_SAVEPATH = []; VFIG_SAVEPATH = []; SORTEDSPIKES_SAVEPATH = [];         DAQSAVEPATH = [];
s_RIGSPECIFIC_SPIKESORT_CONFIG % load PATHinfo

% look up DAQchns
[temp temp temp sexptname DAQchns] = querySpkSortTable(TABLES_SAVEPATH,spkSortFN);
clear temp;

% lookup conditions for each spike
[spikesortdata rawdatainfo] = putStimCondSpkSortData(spkSortFN);


%  VStim Variables during for each DAQfile recording is saved to a table.

if(~exist('cond','var'))|| isempty(cond); cond = struct('notempty',1); end; % will analyze all files that match first STF and meet AnalType specs

[Vstimcellout] = updateVStimCondTable(TABLES_SAVEPATH,{STF.filename});


lastVarParam =[];
for ic = 1:length(cond) % for each condition
    if isstr(AnalType)
        % GET fileind of rawdata that matches cond required
        switch(AnalType)
            case 'ContrastO'
                % *TOD deal with cae of Contrast and fixed orientation
                %         'VarParm1'  'VarParm2'
                indreq = [2 4; 2 4]; typereqcol ={'Orientation','Contrast';'Contrast','Orientation'};
                %            'VarParmVal1'   'VarParmVal2' 'StimulusName'  'StimDuration'    'Baseline'
                indReqByFirstFile = [3 5 6 7 8]; % specified based on first file meeting indreq and conditions
                
            case 'Orientation'
                %         'VarParm1'  'VarParm2'    'VarParmVal2'
                indreq = [2 4 5];        typereqcol ={'Orientation',NaN,[]}; % required for this type of analysis
                %            'VarParmVal1'  'StimulusName'  'StimDuration'    'Baseline'
                indReqByFirstFile = [3 6 7 8]; % specified based on first file meeting indreq and conditions
        end
        % Check that filenames have the same Stim Condition. Below are the
        % parameters that must match to be analyzed together
        [indrawdatafilecond lstrawfilename ConditionDesc indSTFfilecond]= getfiles(cond(ic),indreq,indReqByFirstFile,typereqcol,Vstimcellout,STF,rawdatainfo);
        if ~isempty(indrawdatafilecond)
            [spikeData swcond Nspikes]= concatenateSpiketimesAcrossFiles(spikesortdata,indrawdatafilecond,rawdatainfo,includesweeps);
            totalUNITs = unique(spikeData(:,3)); totalUNITs = totalUNITs(totalUNITs>0);
            
            if ~exist('UNIT','var')|| isempty(UNIT); UNIT =  0; end % default display all units
            if ~all(ismember(UNIT,totalUNITs)); printf('Requested units do not exist in spikesortdatafile'); return ;end
            
            % pre declare for analysis/plotting functions
            maxsample = rawdatainfo(indrawdatafilecond(1)).maxsample;
            dt = rawdatainfo(indrawdatafilecond(1)).dt;
            maxtrigger = rawdatainfo(indrawdatafilecond(1)).maxtrigger;
            Unitcmap = jetm(length(totalUNITs));
            
            if ~isempty(indrawdatafilecond)
                switch(AnalType)
                    case 'ContrastO'
                        ContrastOAnal();
                        
                    case 'Orientation'
                        
                        OrientationAnal();
                        fids = 1000+ic;
                        
                        %                     summarizeTuningofUnits(UNIT,spkSortFN,lstrawfilename,fids);
                        
                end
            end
        else
            printf('** No rawdata found for cond:\n')
            Vstimcellout
        end
    else % this is the original anlysis function. right now it can only be used on the FIRST file in STF
        indrawdatafilecond = 1; indSTFfilecond = indrawdatafilecond;
        [spikeData swcond]= concatenateSpiketimesAcrossFiles(spikesortdata,indrawdatafilecond,rawdatainfo,includesweeps);
        % pre declare for analysis/plotting functions
        maxsample = rawdatainfo(indrawdatafilecond(1)).maxsample;
        dt = rawdatainfo(indrawdatafilecond(1)).dt;
        maxtrigger = rawdatainfo(indrawdatafilecond(1)).maxtrigger;
        [VarParam SParam] = formatVStimCellInfo(Vstimcellout,indSTFfilecond+1); % get the Vstim params of first file
        
        if length(VarParam)==1 && AnalType>0 % catch appropriate cases
            warning('cannot use this AnalType. only 1 VarParam exists')
            AnalType = 0;
        end
        
        %         sAnn = strfind(sf,'\');sAnn = sf(sAnn(end)+1:end);
        nClast = []; TRIG = [];indNaN = [];
        
        [junk TRIG] = getStimCond({fullfile(DAQSAVEPATH,rawdatainfo(indrawdatafilecond(1)).filename)});
        
        FigureName = [VarParam(1).Name SParam.StimName];
        warning('************************** This part not full working yet .... ');
        
        switch AnalType
            case 1 %  plot V1 for each V2
                NVAR = length(VarParam(2).Val);
                for kk = 1:NVAR
                    [junk TRIG] = getStimCond({fullfile(DAQSAVEPATH,rawdatainfo(indrawdatafilecond(1)).filename)});
                    notNAN = ~isnan(TRIG); % have to separate nans because can't do ~ operation on array with NaN values
                    temp = TRIG(notNAN);
                    temp((rem(temp,length(VarParam(1).Val))+ ...
                        ~(rem(temp,length(VarParam(2).Val)))*length(VarParam(2).Val))...
                        ~= kk) = nan; % select conditions for VarParam(2).Val(kk)
                    TRIG(notNAN) = temp; % put back into TRIG
                    TRIG = (floor((TRIG-1)/length(VarParam(2).Val))+1); % rescale conditions to be 1:NV1
                    tempspikeDatacond = spikeData(:,4) ;
                    tempspikeDatacond((rem(tempspikeDatacond,length(VarParam(1).Val))+ ...
                        ~(rem(tempspikeDatacond,length(VarParam(2).Val)))*length(VarParam(2).Val))...
                        ~= kk) = nan; % select conditions for VarParam(2).Val(kk)
                    tempspikeDatacond =  (floor((tempspikeDatacond-1)/length(VarParam(2).Val))+1);
                    
                    Nfig = 100+kk*100;
                    sDesc = sprintf('%s %s', VarParam(2).Name, num2str(VarParam(2).Val(kk)));
                    Vs = sprintf('%s %s %s ', VarParam(2).Name(1:3), num2str(VarParam(2).Val(kk)),VarParam(1).Name(1:3));
                    
                    Val = VarParam(1).Val';
                    helper_CondPLOT
                    slegend{kk} = sDesc;
                end
                legend(slegend);
            case 2 %  plot V2 for each V1
                NVAR = length(VarParam(1).Val);
                for kk = 1:NVAR
                    [junk TRIG] = getStimCond({fullfile(DAQSAVEPATH,rawdatainfo(indrawdatafilecond(1)).filename)});
                    TRIG((~((TRIG<=kk*length(VarParam(2).Val)) & (TRIG>(kk-1)*length(VarParam(2).Val)))))=NaN;
                    Nfig = 100+kk*100;
                    TRIG = TRIG-(kk-1)*length(VarParam(2).Val); % rescale conditions to be 1:NV2
                    tempspikeDatacond =  spikeData(:,4);
                    tempspikeDatacond((~((tempspikeDatacond<=kk*length(VarParam(2).Val)) & (tempspikeDatacond>(kk-1)*length(VarParam(2).Val))))) =  NaN;
                    tempspikeDatacond =  tempspikeDatacond-(kk-1)*length(VarParam(2).Val);
                    
                    
                    sDesc = sprintf('%s %s', VarParam(1).Name, num2str(VarParam(1).Val(kk)));
                    Vs = sprintf('%s %s %s ', VarParam(1).Name(1:3), num2str(VarParam(1).Val(kk)),VarParam(2).Name(1:3));
                    
                    Val = VarParam(2).Val';
                    helper_CondPLOT
                    slegend{kk} = sDesc;
                    
                end
                legend(slegend);
                
                
            case 3 % V1 collapse all V2
                i=1; Nfig = 100+i*100;
                TRIG =  (floor((TRIG-1)/length(VarParam(2).Val))+1);
                tempspikeDatacond= (floor((spikeData(:,4)-1)/length(VarParam(2).Val))+1);
                sDesc = sprintf('%s for all %s', VarParam(1).Name, VarParam(2).Name);
                Vs = VarParam(1).Name(1:3);
                Val = VarParam(1).Val';
                helper_CondPLOT
                
                warning('within parameter condition order may be wrong')
            case 4 %V2 collapse all V1
                i=1; Nfig = 100+i*100;
                TRIG = rem(TRIG-1,length(VarParam(2).Val))+1;
                tempspikeDatacond =rem(spikeData(:,4) -1,length(VarParam(2).Val))+1;
                sDesc = sprintf('%s for all %s', VarParam(2).Name, VarParam(1).Name);
                Vs = VarParam(2).Name(1:3);
                Val = VarParam(2).Val';
                helper_CondPLOT
            otherwise
                Nfig = 1;
                tempspikeDatacond = spikeData(:,4);
                if length(VarParam)==2
                    sDesc = sprintf('%s and %s', VarParam(2).Name, VarParam(1).Name);
                else
                    sDesc = sprintf('%s', VarParam(1).Name);
                end
                
                Vs =[]; for i = 1:length(VarParam);        Vs =  [Vs ' ' VarParam(i).Name(1:3)]; end
                
                if length(VarParam)==2
                    Val = [];
                    
                    for i = 1:VarParam(1).N* VarParam(2).N
                        V1ind = floor((i-1)/VarParam(2).N)+1;
                        V2ind = mod(i,VarParam(2).N)+~mod(i,VarParam(2).N)*VarParam(2).N;
                        Val = [Val; VarParam(1).Val(V1ind) VarParam(2).Val(V2ind)];
                    end
                else
                    Val = VarParam(1).Val';
                end
                
                helper_CondPLOT
        end
        
        
    end
end





%% FIX unit numbesr

% UNIT = 53;


    function ContrastOAnal
        bplot = 1; % for debuggin
        %         % take care of case where VarParm(1) is Contrast no Orientation.
        tagfileind = zeros(size(indrawdatafilecond));
        for i3 = 1:length(indrawdatafilecond)
            [VarParam SParam] = formatVStimCellInfo(Vstimcellout,indSTFfilecond(i3)+1); % get the Vstim params of first file
            if strcmp(VarParam(2).Name,'Contrast');                tagfileind(i3) = 1;            end
        end
        %         % TEMP debug re numbering conditions
        %         temp  =  spikeData(:,4);
        %                 V1 = ceil(temp/length(VarParam(2).Val));
        %                 spikeData(:,4) = ceil(spikeData(:,4)/length(VarParam(2).Val));
        %
        % %                 spikeData(:,4) =  spikeData(:,4)
        % -(V1-1)*length(VarParam(2).Val);
        %  % end of debug
        % make look up table to so that swcond can be switched (in case
        % where Orientation and Contrast are switched in order)
        if all(tagfileind==1) % if there is a switch
            lookupcond = helpswcond([1:length(VarParam(2).Val)*length(VarParam(1).Val)],length(VarParam(2).Val),length(VarParam(1).Val));
            for i3 = 1:length(indrawdatafilecond)
                cumNspikes = []; % switch conditions for files used in analysis (indrawdatafilecond)
                if isempty(cumNspikes);    for i4 = 1: length(Nspikes); cumNspikes = [0 cumsum(Nspikes)]+1; end ;end
                temp = spikeData(cumNspikes(i3):cumNspikes(i3+1)-1,4);
                spikeData(cumNspikes(i3):cumNspikes(i3+1)-1,4) = helpswcond(temp,length(VarParam(2).Val),length(VarParam(1).Val));
            end
            % flip naming
            temp = VarParam(1).Name;     VarParam(1).Name = VarParam(2).Name; VarParam(2).Name = temp;
            temp = VarParam(1).Val;     VarParam(1).Val = VarParam(2).Val; VarParam(2).Val = temp;
        elseif ~(all(tagfileind==1)||all(tagfileind==0)); error('Cannot analyze different files that are Var1= Contrast,  Var2 = Orientation, AND Var1 = Orientation , Var2 = Contrast');
        else lookupcond = [1:length(VarParam(2).Val)*length(VarParam(1).Val)];
        end
        
        
        if ~isempty(lastVarParam)
            % new figure for Var Param is different condition
            if ~isequal(VarParam,lastVarParam);                    fid = 10000+100*(ic-1);                 end
        else                fid = 20000;            end
        
        lastVarParam = VarParam;
        
        FigureName = ['Contrast: ' SParam.StimName];
        
        % define nCond (after they are switched)
        nCond(1) = length(VarParam(1).Val);
        nCond(2) = length(VarParam(2).Val);
        CondVal = VarParam(2).Val; % orientations
        CondVal2 = VarParam(1).Val; % contrast
        
        
        PSTHcond = nan(length(UNIT),nCond(1), nCond(2),maxsample-2);
        %          smoothPSTHcond(jj,iVar1,iVar2,:)
        smoothPSTHcond = nan(length(UNIT),nCond(1), nCond(2),maxsample-2);
        PSTHnSw = nan(length(UNIT),nCond(1),nCond(2)); % number of sweeps/trials in psth condition
        xtime = [1:size(PSTHcond,4)]*dt;
        WOI = [0 maxsample-3];
        % define epocs to analyze
        if isempty(analyzewindow) % use user specified window if not empty
            sampWIND = round([0 min(size(PSTHcond,4)-1,SParam.StimDuration/dt)]+1);
            analyzewindow = round([0 min((size(PSTHcond,4)-1)*dt,SParam.StimDuration)]);
        else
            sampWIND = round([analyzewindow(1)/dt min(size(PSTHcond,4)-1,analyzewindow(2)/dt)]+1);
        end
        % baseline periods is when gray background is on after
        % 100ms of the stimulus being off
        if ~isempty(SParam.Baseline) && SParam.Baseline
            baselinesampWIND = round( [(SParam.StimDuration+0.1)/dt min(size(PSTHcond,3),(SParam.StimDuration+SParam.Baseline)/dt)]);
        else baselinesampWIND = [NaN NaN]; end
        
        % baseline periods is when gray background is on after
        % 100ms of the stimulus being off
        if ~isempty(SParam.Baseline) && SParam.Baseline
            baselinesampWIND = round( [(SParam.StimDuration+0.1)/dt min(size(PSTHcond,3),(SParam.StimDuration+SParam.Baseline)/dt)]);
        else baselinesampWIND = [NaN NaN]; end
        % setup figure subplots
        nCol = max(nCond(2),4);
        if mod(nCol,2); nCol = nCol+1;  end  % make even makes thinks easier
        totsubplot = nCol * nCond(1) + 2*nCond(1);  % totla number of subplots
        nRow = ceil(totsubplot/nCol);
        % PLOT Positions BA
        nRowrast  = nCond(1);
        nColrast  =  nCond(2);
        axparams.matpos = [0 0.02 0.75 0.75];%         offset[left top]) , fraction [width height] of remaining figure (i.e. 1-(FIGMARG)  that is taken up by matrix(measured from WITHIN margin of figure)
        axparams.figmargin = [0 0 0.010 0.01]; % margins within whole figure [ LEFT RIGHT TOP BOTTOM]
        axparams.matmargin = [.03 0 0.0 0]; % margins within? whole matrix [ LEFT RIGHT TOP BOTTOM] note MATMARG, can also be thought of as position
        axparams.cellmargin = [0.01 0.015 0.01 0.01]; % MARGIN around each axies (defines offset and size of each axes with a cell)
 
        % bottom plots of CRF and ori
        axparamsbot.matpos = [ 0 0.75 1 0.25];
        axparamsbot.figmargin =  axparams.figmargin;
        axparamsbot.matmargin = [.03 0.03 0.04 0.03];
        axparamsbot.cellmargin = [0.01 0.12 0.00 0.00];
        
        % spk plots 
        axparamsSpk.matpos = [.75 0.15 .25 0.5];
        axparamsSpk.figmargin =  axparams.figmargin;
        axparamsSpk.matmargin = [.03 0.05 0.01 0.08];
        axparamsSpk.cellmargin = [0.01 0.01 0.05 0.05];

        for jj=1:length(UNIT) %each unit
            
            ifid = fid+jj*10;
            figure(ifid);  orient tall;      set(ifid,'position',[2          39        1178         919]); % MAKE rig specific
            % choose raster or psth
            if 1  % Rasterplot
                clear h ph
                for iVar1 = 1:nRowrast % ROW 
                    for iVar2 = 1:nColrast % COL
                        icond = (iVar1-1)*nCond(2)+iVar2;
                        maskout = maskspikes(UNIT(jj),icond,spikeData);
                        h(icond) = axesmatrix(nRowrast,nColrast,(iVar1-1)*nColrast + iVar2,axparams);
                        set(h(icond),'XAxisLocation','bottom',...
                            'YAxisLocation','left',...
                            'Color','none',...
                            'XColor','k','YColor','k');
                        
                        spiketime = (spikeData(maskout,1)+(spikeData(maskout,2)-1)*maxsample); %convert to continous spiketimes across sweeps
                        if iVar1==1;                            stemp = '';  if iVar2==1;stemp = 'Orient';end;  stemp = sprintf('%s %s',stemp,num2str(CondVal(iVar2))  );title(stemp); end
                        if iVar2==1;                   ylabel(num2str(   CondVal2(iVar1),'%1.2f'));                 end
                        plotset(1);set(gca,'box','on');hold all;
                        ph{icond} = rasterplot(spiketime,dt/1e3,maxsample,1,ifid,.95,[]);axis tight;
                    end
                    
                end
                if ~isempty(h); %                     format axis % auto increment color
                    h = h(h>0); % remove empty plots that didn't get a handle
%                     stemp = sprintf('\\bf{Unit} %s\ntrial',num2str(UNIT(jj)));
                    stemp = sprintf('Contrast\n%s',get(get(h(1),'YLabel'),'String'));
                     setColorbyColorOrder(h,ph);
                    set(get(h(1),'YLabel'),'String',stemp);
                    
                    ind = nColrast*(nRowrast-1)+1; % index of bottom left axes (which will have x axis info)
                    
                    set(h,'YTickLabel',[],'YTick',[]);
                    set(h(h~=h(ind)),'XTickLabel',[]);
                    set(get(h(ind),'XLabel'),'String','sec');
                    temp  = get(get(h(ind),'XLabel'),'Position');
%                     set(get(h(ind),'XLabel'),'Position',temp.*[1 0.75 1]); % bring the xlabel closer to axis
                    temp = ([min(cellfun(@min,get(h,'YLIM'))) max(cellfun(@max,get(h,'YLIM')))]);
                    linkprop(h,{'box','Color','YLIM','XLIM'});
                    set(h,'YLIM',temp,'XLIM',[ 0 maxsample*dt],'box','off');

                    % analysis period
                    set(gcf,'CurrentAxes',h(1));
                    line ([sampWIND(1)*dt sampWIND(2)*dt], [1 1].*max(YLIM)*1.03,'linewidth',2,'color','r')
                    % stimulus period
                    line ([0 SParam.StimDuration], [1 1].*max(YLIM)*1.06,'linewidth',2,'color','k'); set(gca,'YLIM',[min(ylim) 1.1*max(ylim)]);
                    
                    
                end
            end
            if 1  % % PSTH
                for iVar1 = 1:nCond(1) % contrast
                    for iVar2 = 1:nCond(2) % orietation
                        
                        %%% calculate
                        icond = (iVar1-1)*nCond(2)+iVar2;
                        maskout = maskspikes(UNIT(jj),icond,spikeData);
                        eventtime = (find(swcond==find(lookupcond==icond))-1)*maxsample+2; % eventtime is the beginning of each sweep (with a few sample offset), (must convert swcond number if Var1 and Var2 are switched)
                        spiketime = (spikeData(maskout,1)+(spikeData(maskout,2)-1)*maxsample); %convert to continous spiketimes across sweeps
                        [eventHist skipped] = perEventSpkHist(eventtime, spiketime, maxsample, WOI);
                        if ~isempty(eventHist)
                            if size(eventHist,1)>1
                                PSTHcond(jj,iVar1,iVar2,:) = sum(eventHist);
                            else
                                PSTHcond(jj,iVar1,iVar2,:) = (eventHist);
                            end
                            if(any(skipped)); error('should not be skipping'); end % number of sweeps for each unit for each condition that went into PSTH
                            
                            PSTHnSw(jj,iVar1,iVar2,:) = length(eventtime);
                        end
                        %%% plot
                        if bplot
                            % TODO 1) TOO slow!! 2)replace with convoltion? add bootstrap confidence
                            if sum(PSTHcond(jj,iVar1,iVar2,:))>0 % save time
                                %                         tic
                                smoothPSTHcond(jj,iVar1,iVar2,:) = filterdata(single(squeeze(PSTHcond(jj,iVar1,iVar2,:))),dt,5,0);%          b = smooth(a,.01/dt);
                                %                         toc
                            else                             smoothPSTHcond(jj,iVar1,iVar2,:) = nan(1,1,size(PSTHcond(jj,iVar1,iVar2,:),4)); end
                            
                             hp(icond) = axes('Position',get(findobj(h(icond),'Type','axes'),'Position'),...
                                'XAxisLocation','top',...
                                'YAxisLocation','right',...
                                'Color','none',...
                                'XColor','k','YColor','r');
                            if ~all(isnan(smoothPSTHcond(jj,iVar1,iVar2,:)))
                                line(xtime,squeeze(smoothPSTHcond(jj,iVar1,iVar2,:))/dt/PSTHnSw(jj,iVar1,iVar2,:),'linewidth',2,'parent',hp(icond),'color','r','linestyle','-');
                            else line(nan,nan,'parent',hp(icond));end
                            hold all; plotset(1);
                            axis tight;
                            % baseline period baselinesampWIND
                            line(baselinesampWIND(1)*dt.*[1 1],[0 5],'color',[0 0 0],'parent',hp(icond))
                        end
                    end
                end
                
                if ~isempty(hp); %                     format axis
                    if bplot
                        hp = hp(hp>0); % remove empty plots that didn't get a handle
                        for i = 1:length(hp); set(hp(i),'Position',get(findobj(h(i),'Type','axes'),'Position')); end
                        if length(hp)>1
                            temp = ([min(cellfun(@min,get(hp,'YLIM'))) max(cellfun(@max,get(hp,'YLIM')))]);
                            linkprop(hp,{'box','Color','YLIM'});
                            linkprop([h hp],'XLIM');
                            set(hp,'YLIM',temp,'XLIM',[ 0 maxsample*dt],'box','off');
                        end
                        ind = nColrast;% index of axes on far right of first row
                        set(hp(1:end),'XTickLabel',[],'XTick',[]);
                        set(hp(hp~=hp(ind)),'YTickLabel',[]);
                        set(get(hp(ind),'YLabel'),'String','Rate(Hz)')
                        
                    end
                end
            end
            % PLOT **CRF** contrast response functino for each orientation
            axesmatrix(1,2,1,axparamsbot);
            cmap = repmat(linspace(0.85,0,nCond(2))',1,3); % grayscale color map
            sl2 = cell(iVar2,1); clear sp;
            for iVar2 = 1:nCond(2)
                spikeRate = sum(PSTHcond(jj,:,iVar2,sampWIND(1):sampWIND(2)),4)./PSTHnSw(jj,:,iVar2,:)/(diff(sampWIND)*dt);
                
                sp(iVar2) =  plot(CondVal2,spikeRate,'-o','color',cmap(iVar2,:),'linewidth',2); hold on;
                spikeRateBL = sum(PSTHcond(jj,:,iVar2,baselinesampWIND(1):baselinesampWIND(2)),4)./PSTHnSw(jj,:,iVar2,:)/(diff(baselinesampWIND)*dt);
                
                plot(CondVal2,spikeRateBL,'--','color',cmap(iVar2,:),'linewidth',2); hold on;
                sl2{iVar2} = num2str(CondVal(iVar2));
            end
            
            % find condition with maximum number of spikes and fit it with
            % c50
            sPSTH= squeeze(sum(PSTHcond(jj,:,:,sampWIND(1):sampWIND(2)),4)./PSTHnSw(jj,:,:,:)/(diff(sampWIND)*dt))';
            [temp indexmax] = max(sum(sPSTH'));
            
            %             mcrf = mean(sum(PSTHcond(jj,:,:,sampWIND(1):sampWIND(2)),4)./PSTHnSw(jj,:,:,:)/(diff(sampWIND)*dt),3);
            %             plot(CondVal2,sPSTCH(:,h),'-r','linewidth',2);
            % fit crf
            if iVar2>1; temp = sPSTH(indexmax,:)'; else temp = sPSTH; end % special case for only 1 orientation
            [c2 gof2] = fitcontrast(CondVal2',temp);
            crf_c50 = c2.c50;  crf_n = c2.n; crf_rsq = gof2.rsquare;
            plot(c2,'b') ;
            stemp = sprintf('c50: %1.1f, n: %1.1f, r^2: %1.2f',crf_c50,crf_n,crf_rsq);
            title(stemp)
            axis tight; hl = legend(sp,sl2,'box','off','color','none','location',[0.4042    0.08    0.07    0.15]);           xlabel('Contrast'); ylabel('rate (hz)');            plotset(1);
            
            %
            % PLOT **orientation tuning** for each contrast
            if iVar2>1; % only if more than 1 orientation
                h =  axesmatrix(1,2,2,axparamsbot);
                cmap = repmat(linspace(0.85,0,nCond(1))',1,3); % grayscale color map
                sl2 = cell(iVar1,1);
                for iVar1 = 1:nCond(1)
                    spikeRate = squeeze(sum(PSTHcond(jj,iVar1,:,sampWIND(1):sampWIND(2)),4)./PSTHnSw(jj,iVar1,:,:)/(diff(sampWIND)*dt));
                    
                    plot(CondVal,spikeRate,'-o','color',cmap(iVar1,:),'linewidth',2); hold on;
                    sl2{iVar1} = num2str(CondVal2(iVar1),'%1.2f');
                end
                axis tight;       xlabel('Orient');      set(h,'XTick',sort(CondVal)); plotset(1); hl = legend(sl2,'box','off','color','none','location',[0.8765    0.08    0.07    0.17]);
            end
            % plot waveform
             hp(1) =  axesmatrix(2,1,1,axparamsSpk);
             hp(2) =  axesmatrix(2,1,2,axparamsSpk);
            [waves Nspks NBAD] = spikequality(UNIT(jj),spkSortFN,hp);
            t = dt*[1:size(waves,2)]*1e3; % ms
            
            sannotate = helpAnnotate(jj,ifid,FigureName,SParam);  %NESTED function (shares same workspace)
            
            if bsave
                % save figure
                diranaltype = 'Crf';
                sexptname_score = regexprep(sexptname,'\s','_');
                figfname = sprintf('crf_%sC%s_U%d_INT%s_%s%s%s', getfileindex(STF),num2str(cell2mat(DAQchns)),UNIT(jj),num2str(analyzewindow(1)*10,'%d'),num2str(analyzewindow(2)*10,'%d'),regexprep(SParam.StimName,'\s','_'),regexprep(SParam.sDesc,'\s','_'));
                
                D = fullfile(VFIG_SAVEPATH,diranaltype);
                if isempty(dir(fullfile(D,sexptname_score)));mkdir(D,sexptname_score); end
                figdir =fullfile(D,sexptname_score);
                savefigure(figdir,[],figfname,'fig',ifid)
                
                % Save parameters
                unitsummary.ExptName = sexptname;
                unitsummary.rawdata_in_analysis = {rawdatainfo(indrawdatafilecond).filename};
                unitsummary.DAQchns = DAQchns;
                unitsummary.Unit = UNIT(jj);
                unitsummary.spkSortFN = spkSortFN;
                unitsummary.crf_c50 = crf_c50;
                unitsummary.crf_n = crf_n;
                unitsummary.crf_rsq = crf_rsq;
                unitsummary.crffit_c = c2;
                unitsummary.crffit_gof = gof2;
                unitsummary.spf_Freq = SParam.spfFreq ;
                unitsummary.PSTH = squeeze(PSTHcond(jj,:,:,:));
                unitsummary.PSTHnSw = squeeze(PSTHnSw(jj,:,:,:));
                unitsummary.PATH_ContrastO_figure = fullfile(diranaltype,sexptname_score,figfname);
                
                unitsummaryCond.label = label; % label designed to be helpful in identifying this
                unitsummaryCond.maxsample = maxsample;
                unitsummaryCond.maxtrigger = maxtrigger;
                unitsummaryCond.dt = dt;
                unitsummaryCond.CondVal = CondVal;
                unitsummaryCond.CondVal2 = CondVal2;
                unitsummaryCond.sampWIND = sampWIND;
                unitsummaryCond.baselinesampWIND = baselinesampWIND;
                unitsummaryCond.totalUNITs = totalUNITs;
                unitsummaryCond.sannotate = sannotate;
                unitsummaryCond.ConditionDesc = ConditionDesc; % important , defines conditions that rawdatafile met to be included in analysis
                unitsummaryCond.includesweeps = includesweeps; % important , defines sweeps in rawdatafile  included in analysis
                
                %save to file
                savetuningProperties(unitsummary,unitsummaryCond)
            end
        end
    end

    function OrientationAnal
        [VarParam SParam] = formatVStimCellInfo(Vstimcellout,indSTFfilecond(1)+1); % get the Vstim params of first file
        % TO DO ** should not display this though as it could be
        % missleading (if there are multiple files that wouldn't
        % necessary have all the same SParam , but just those specifed
        % by cond
        %                  fid = 100;
        fid = 100+100*(ic-1);
        %         if ~isempty(lastVarParam)
        %             % new figure for Var Param is different condition
        %             if ~isequal(VarParam,lastVarParam);                    fid = 100+100*(ic-1);                 end
        %         end
        %         lastVarParam = VarParam;
        %
        FigureName = ['Ori: ' SParam.StimName];
        nCond = length(VarParam(1).Val);
        CondVal = VarParam(1).Val;
        cmap = jetm(nCond);
        
        % predefine variables
        PSTHcond = nan(length(UNIT),nCond,maxsample-2);
        smoothPSTHcond = nan(length(UNIT),nCond,maxsample-2,1);
        PSTHnSw = nan(length(UNIT),nCond,1); % number of sweeps/trials in psth condition
        xtime = [1:size(PSTHcond,3)]*dt;
        WOI = [0 maxsample-3];
        
        % define epocs to analyze
        sampWIND = round([0 min(size(PSTHcond,3)-1,SParam.StimDuration/dt)]+1);
        % baseline periods is when gray background is on after
        % 100ms of the stimulus being off
        if ~isempty(SParam.Baseline) && SParam.Baseline
            baselinesampWIND = round( [(SParam.StimDuration+0.1)/dt min(size(PSTHcond,3),(SParam.StimDuration+SParam.Baseline)/dt)]);
        else baselinesampWIND = [NaN NaN]; end
        
        % setup figure subplots
        nCol = max(nCond,4);
        if mod(nCol,2); nCol = nCol+1;  end  % make even makes thinks easier
        totsubplot = nCol * 2 + nCol*2;  % totla number of subplots
        nRow = ceil(totsubplot/nCol);
        %         clear PSTHcond PSTHnSw smoothPSTHcond
        for i=1:length(UNIT) %each unit
            ifid = fid+i;
            figure(ifid); orient landscape
            
            
            % Rasterplot
            helpRasterplot(i,ifid,nRow,nCol,nCond,CondVal); %NESTED function (shares same workspace)
            
            %%% PSTH ( hard to make into function because PSTHcond need to
            %%% be passsed back and for which are big and can get slo
            clear h
            for j = 1:nCond
                %%% calculate
                maskout = maskspikes(UNIT(i),j,spikeData);
                eventtime = (find(swcond==j)-1)*maxsample+2; % eventtime is the beginning of each sweep (with a few sample offset)
                spiketime = (spikeData(maskout,1)+(spikeData(maskout,2)-1)*maxsample); %convert to continous spiketimes across sweeps
                [eventHist skipped] = perEventSpkHist(eventtime, spiketime, maxsample, WOI);
                if ~isempty(eventHist)
                    if size(eventHist,1)>1
                        PSTHcond(i,j,:) = sum(eventHist);
                    else
                        PSTHcond(i,j,:) = (eventHist);
                    end
                    if(any(skipped)); error('should not be skipping'); end % number of sweeps for each unit for each condition that went into PSTH
                    
                    PSTHnSw(i,j,:) = length(eventtime);
                end
                % TODO replace with convoltion? add bootstrap
                % confidence
                if sum(PSTHcond(i,j,:))>0 % save time
                    smoothPSTHcond(i,j,:) = filterdata(single(squeeze(PSTHcond(i,j,:))),dt,5,0);%          b = smooth(a,.01/dt);
                else                             smoothPSTHcond(i,j,:) = nan(1,1,size(PSTHcond(i,j,:),3)); end
                
                
                %%% plot
                h(j) = subplot(nRow,nCol,j+nCol);
                lh = line(xtime,squeeze(smoothPSTHcond(i,j,:))/dt/PSTHnSw(i,j,:),'linewidth',2,'tag','data'); hold all; plotset(1);
                
                Co = get(h(j) ,'colorOrder'); Cind = length(findobj(get(gca,'Children'),{'Type','Tag'},{'line','data'}));
                helpsetcolor({lh},Co(mod(Cind,size(Co,1))+size(Co,1)*~mod(Cind,size(Co,1)),:))
                axis tight;
                % baseline period baselinesampWIND
                line(baselinesampWIND(1)*dt.*[1 1],[0 5],'color',[0 0 0])
            end
            if ~isempty(h); %                     format axis
                stemp = sprintf('\\bf{Unit} %s\n spk rate (Hz) **CHECK',num2str(UNIT(i)));
                set(get(h(1),'YLabel'),'String',stemp,'interpreter','Latex','color',Unitcmap(UNIT(i),:));
                set(get(h(1),'XLabel'),'String','sec');
                set(h(2:end),'XTickLabel',[],'YTickLabel',[]);
                if length(h)>1
                temp = ([min(cellfun(@min,get(h,'YLIM'))) max(cellfun(@max,get(h,'YLIM')))]);
                linkprop(h,{'box','Color','YLIM','XLIM'});
                xlim([ 0 maxsample*dt]);ylim(temp);
                end
            end
            
            % Summarize across tuning properties
            
            % calculate  tuning
            spikeRate = sum(PSTHcond(i,:,sampWIND(1):sampWIND(2)),3)./PSTHnSw(i,:,1)/(diff(sampWIND)*dt);
            r = [spikeRate'; spikeRate(1)]; r(isnan(r)) = 0; %replace nan with zero WATCH out this could be bad..
            theta = [CondVal'; CondVal(1)];
            if ~isnan(baselinesampWIND) % baseline calculated for all conditions
                baselinespikeRate = sum(PSTHcond(i,:,baselinesampWIND(1):baselinesampWIND(2)),3)./PSTHnSw(i,:,1)/(diff(baselinesampWIND)*dt);
                rBL = [baselinespikeRate'; baselinespikeRate(1)]; rBL(isnan(rBL)) = 0;
            else                       rBL = nan(size(theta));                    end
            
            
            % calculate tuning metrics
            [tuntemp indmax]= calcDirnOriTuning(r,theta);
            prefDirn = theta(indmax(1));
            % get F1,F0 Response
            F1F0 = getF1F0(squeeze(PSTHcond(i,indmax(1),sampWIND(1):sampWIND(2))),dt,SParam.tempFreq );
            % get max and median Stimulus/Baseline rate
            temp = squeeze(smoothPSTHcond(i,indmax(1),sampWIND(1):sampWIND(2)))./PSTHnSw(i,indmax(1),1)/dt;
            rateMax(1) = max(temp(:));
            rateMed(1) = median(temp(:));
            if ~isnan(baselinesampWIND)
                temp = nansum(smoothPSTHcond(i,:,baselinesampWIND(1):baselinesampWIND(2)),2)/sum(PSTHnSw(i,:,1))/dt;
                rateMax(2)= max(temp(:));
                rateMed(2)=  median(temp(:));
            else rateMax(2)= nan;  rateMed(2) = nan; end
            
            % plotting tunign properties
            clear htemp
            subplot(nRow,nCol,nCol*2.*[1 1]+[1 nCol/2]); % subplot takes half the row % polar plot
            htemp(1) = plot(theta(1:end-1),r(1:end-1),'-o'); hold all; plotset(1)
            set(htemp,'linewidth',2);
            htemp(2) = plot(theta(1:end-1),rBL(1:end-1),'-');  plotset(1) % baseline
            setColorbyColorOrder(gca,{htemp},2)
            ylabel('rate (Hz)');
            set(gca,'XTick',[0:45:360])
            
            subplot(nRow,nCol,nCol*2.*[1 1]+[nCol/2+1 nCol]); % subplot takes the next half
            htemp = polar(deg2rad(theta),r/max(r)); hold all
            set(htemp,'linewidth',2)
            
            clear h htemp;
            h = subplot(nRow,nCol,nCol*3.*[1 1]+[nCol/2+1 nCol]);
            htemp(1) = plot(ones(size(tuntemp)),tuntemp.OSI,'o');hold all;
            % Dirn selectivity                                            % (Rpref-Roppsite)/(Rpref+Roppsite)
            htemp(2)= plot(2*ones(size(tuntemp)),tuntemp.DSI,'o');
            htemp(3)= plot(3*ones(size(tuntemp)),tuntemp.VAVG,'o'); % vector average
            htemp(4)= plot(4*ones(size(F1F0)),F1F0,'o'); % vector average
            set(gca,'XTick',[1:4],'XTickLabel',{'Osi','Dsi','Os_Vec','F1F0'}); xlim([0.5 4+.5]);
            temp = ylim; ylim([0 1])
            ylabel('$$\frac{R_p- R_o}{R_p+ R_o}$$','Interpreter','latex');   plotset(1)
            setColorbyColorOrder(h,{htemp},4)
            
            clear h htemp; % SOEMTHING WRONG WITH COLOR BUT DONT KNOW WHAT
            h =subplot(nRow,nCol,nCol*3.*[1 1]+[1 nCol/2]);
            indexX = 1;
            htemp(1)= plot(indexX,rateMax(1),'o'); hold all; % stimulus
            htemp(2)= plot(indexX,rateMax(2),'o');% baseline
            indexX = 2;
            htemp(3)= plot(indexX,rateMed(1),'o'); % stimulus
            htemp(4)= plot(indexX,rateMed(2),'o');% baseline
            set(gca,'XTick',[1:indexX],'XTickLabel',{'Max','Med' }); xlim([0.5 indexX+.5]);
            ylabel('rate(Hz)');   plotset(1)
            
            setColorbyColorOrder(h,{htemp},4)
            
            % plot waveform
            [waves Nspks NBAD] = spikequality(UNIT(i),spkSortFN)  ;
            
            
            sannotate = helpAnnotate(i,ifid,FigureName,SParam);  %NESTED function (shares same workspace)
            
            if bsave
                % save figure
                diranaltype = 'Orient';
                sexptname_score = regexprep(sexptname,'\s','_');
                figfname = sprintf('Ori_%sC%s_U%d_%s%s', getfileindex(STF),num2str(cell2mat(DAQchns)),UNIT(i),regexprep(SParam.StimName,'\s','_'),regexprep(SParam.sDesc,'\s','_'));
                D = fullfile(VFIG_SAVEPATH,diranaltype);
                if isempty(dir(fullfile(D,sexptname_score)));mkdir(D,sexptname_score); end
                figdir =fullfile(D,sexptname_score);
                
                savefigure(figdir,[],figfname,'fig',ifid)
                
                unitsummary.ExptName = sexptname;
                unitsummary.rawdata_in_analysis = {rawdatainfo(indrawdatafilecond).filename};
                unitsummary.DAQchns = DAQchns;
                unitsummary.Unit = UNIT(i);
                unitsummary.spkSortFN = spkSortFN;
                unitsummary.max_rate_hz = rateMax(1);
                unitsummary.BL_max_rate = rateMax(2);
                unitsummary.med_rate = rateMed(1);
                unitsummary.BL_med_rate = rateMed(2);
                unitsummary.pref_Dirn = prefDirn;
                unitsummary.spf_Freq = SParam.spfFreq ;
                unitsummary.OSI = tuntemp.OSI;
                unitsummary.DSI = tuntemp.DSI;
                unitsummary.Vavg = tuntemp.VAVG;
                unitsummary.F1F0 = F1F0;
                unitsummary.PSTH = squeeze(PSTHcond(i,:,:));
                unitsummary.PSTHnSw = squeeze(PSTHnSw(i,:,:));
                unitsummary.PATH_Orientation_figure = fullfile(diranaltype,sexptname_score,figfname);
                
                unitsummaryCond.bMask = SParam.bMask;
                unitsummaryCond.maxsample = maxsample;
                unitsummaryCond.maxtrigger = maxtrigger;
                unitsummaryCond.dt = dt;
                unitsummaryCond.CondVal = CondVal;
                unitsummaryCond.sampWIND = sampWIND;
                unitsummaryCond.baselinesampWIND = baselinesampWIND;
                unitsummaryCond.totalUNITs = totalUNITs;
                unitsummaryCond.sannotate = sannotate;
                unitsummaryCond.ConditionDesc = ConditionDesc; % important , defines conditions that rawdatafile met to be included in analysis
                
                savetuningProperties(unitsummary,unitsummaryCond)             %                     save tunig properties
                
            end
            
        end % unit
    end
    function [h ph] = helpRasterplot(i,ifid,nRow,nCol,nCond,CondVal)
        clear h ph;
        for j = 1:nCond
            maskout = maskspikes(UNIT(i),j,spikeData);
            h(j) = subplot(nRow,nCol,j);
            spiketime = (spikeData(maskout,1)+(spikeData(maskout,2)-1)*maxsample); %convert to continous spiketimes across sweeps
            title([num2str(CondVal(j))]); plotset(1);set(gca,'box','on');hold all;
            temp = get(h(j),'Tag'); 
            Cind = str2num(temp(regexp(temp,'\d'))) ; % I use the tag to count the number of raster on the same plot for coloring each one different
            if isempty(Cind); Cind = 1; else; Cind= Cind+1;  end; set(h(j),'Tag',['Raster ' num2str(Cind)]);
            Co = get(h(j) ,'colorOrder');
            temp = Co(mod(Cind,size(Co,1))+size(Co,1)*~mod(Cind,size(Co,1)),:);
            ph{j} = rasterplot(spiketime,dt/1e3,maxsample,1,ifid,.95,temp);axis tight;
        end
        if ~isempty(h); %                     format axis % auto increment color
            stemp = sprintf('\\bf{Unit} %s\ntrial',num2str(UNIT(i)));
            set(get(h(1),'YLabel'),'String',stemp,'interpreter','Latex','color',Unitcmap(UNIT(i),:));
            set(get(h(1),'XLabel'),'String','sec');
            set(h(2:end),'XTickLabel',[],'YTickLabel',[]);
            if length(h)>1
                temp = ([min(cellfun(@min,get(h,'YLIM'))) max(cellfun(@max,get(h,'YLIM')))]);
                linkprop(h,{'box','Color','YLIM','XLIM'});
                xlim([ 0 maxsample*dt]);ylim(temp);
            end
        end
    end
    function sannotate = helpAnnotate(i,ifid,FigureName,SParam)
        % annotate plot % with Vstim conditions, and rawdata file
        stemp = sprintf('  %s U: %d %s',regexprep(sexptname,'\s',' '),UNIT(i),FigureName);
        set(ifid,'Name', stemp);
        sannotate = sprintf('%s\n%s\n%s',regexprep(sexptname,'\s',' '),SParam.StimName,SParam.sDesc);
        plotAnn( sprintf('\n\n\n\n\n%s',sannotate),ifid,1);
        tempcell = {rawdatainfo(indrawdatafilecond).filename}; stemp = [];
        for i2 = 1:length(tempcell); % get indexes of rawdata files in this analysis
            tempind = strfind(tempcell{i2},'_');
            stemp = [stemp tempcell{i2}(tempind(end):end) 'f'];
        end
        sannotate2 = sprintf('files: %s\n%s',stemp,label);                h = plotAnn( sannotate2,ifid,3);
        
        temph = findobj(gcf,'-regexp','Tag','Raster\s\d'); % find any plots with Raster tag and set the text to the same color as sthe raster
        if ~isempty(temph)
            temp = get(temph(1),'Tag');             Cind = str2num(temp(regexp(temp,'\d'))) ; % I use the tag to count the number of raster on the same plot for coloring each one different
            Co = get(temph(1),'colorOrder'); sC = Co(mod(Cind,size(Co,1))+size(Co,1)*~mod(Cind,size(Co,1)),:); % make sure index is not larger then colororder matrix
        else sC = [0 0 1]; end
            
            set(h,'color',sC);
    end

%         ph can be 1xN cell each element with 0 or many handles all handeles are set to have color property C

    function helper_CondPLOT
        C =  unique(TRIG);
        C = C(~isnan(C));
        nC=length(C);
        
        PSTHcond =     nan(length(UNIT),nC,maxsample-2);
        PSTHnSw = nan(length(UNIT),nC); % number of sweeps/trials in psth condition
        
        if nC > length(Val)
            warning('Detected more conditions then specified Values');
            nC =  length(Val);
        end
        F = 0;minsw = 0; % for compatiblity with below version
        
        
        clear Cond;
        for icond=1:nC; Cond(icond).I_sw = find(TRIG==C(icond))+minsw; end % get sweeps of each condition
        
        %%% RASTERPLOT
        if 1
            cmap = jetm(length(UNIT));
            for j = 1:length(UNIT) % make struct for CHRONUX psth
                for icond=1:nC
                    tmp = 0; % number of spikes per condition (per unit)
                    unit(j).Cond(icond).times = [];
                    maskout = maskspikes(UNIT(j),icond,[spikeData(:,[1:3]) tempspikeDatacond]);
                    sum(maskout)
                    temp = (spikeData(maskout,1)+(spikeData(maskout,2)-1)*maxsample); %convert to continous spiketimes across sweeps
                    unit(j).Cond(icond).times  = temp;
                    unit(j).Cond(icond).nspikes = length(temp);
                    
                end
            end;
        end
        
        if 0
            fid = Nfig+F; figure(fid); set(fid,'Name', ['RASTER ' sDesc]); clf; r = 2; clear h;
            tic
            for icond=1:nC
                for j = 1:length(totalUNITs)
                    h(icond) = subplot(r,ceil(nC/r),icond);
                    %%TODO change to manual update plot in rasterplot
                    rasterplot(unit(j).Cond(icond).times/dt/1000,dt,maxsample,1,fid,.95,cmap(j,:));
                    %         rasterplot(UNIT(j).Cond(icond).times,dt,maxsample,1,fid,2*4,cmap(j,:));
                    if icond> 1; axis off; else ylabel('trial'); xlabel('sec'); end
                end
                %              title(['Cond: ' num2str(icond)]); plotset(1);set(gca,'box','on')
                title([Vs ' ' num2str(Val(icond,:))]); plotset(1);set(gca,'box','on')
            end
            toc
            if ~isempty(h);        linkaxes(h);    linkprop(h,{'box','TickDir','Color'}); end
            axis tight; xlim([ 0 maxsample*dt])
            sannotate = helpAnnotate(1,fid,FigureName,SParam);  %NESTED function (shares same workspace)
            
            
            
        end
        %         if bAnalyzeAllFilesTogether;
        %         end
        
        
        %%% PSTH
        % psth using my code
        try
            WOI = [0 maxsample-3];
            for icond=1:nC
                %                 eventtime = ((Cond(icond).I_sw-1)*maxsample+2);
                eventtime = (find(swcond==icond)-1)*maxsample+2; % eventtime is the beginning of each sweep (with a few sample offset)
                for j = 1:length(UNIT)
                    maskout = maskspikes(UNIT(j),icond,[spikeData(:,[1:3]) tempspikeDatacond]);
                    spiketime = (spikeData(maskout,1)+(spikeData(maskout,2)-1)*maxsample)';
                    [eventHist skipped] = perEventSpkHist(eventtime, spiketime, maxsample, WOI);
                    if size(eventHist,1)>1
                        PSTHcond(j,icond,:) = sum(eventHist);
                    else
                        PSTHcond(j,icond,:) = (eventHist);
                    end
                    PSTHnSw(j,icond) = length(eventtime); % number of sweeps for each unit for each condition that went into PSTH
                end
            end
        catch
            j
        end
        
        %plot psth ( in HZ)
        cmap = jetm(nC);
        Unitcmap = jetm(length(UNIT));
        fid = Nfig+10+F;figure(fid);clf;  set(fid,'Name', ['PSTH ' sDesc]); clf; r = 2; clear h;
        %         set(gcf,'Position',[2133         -68         548         872])
        set(gcf,'Position',[1690+Nfig/2        -138         611         466])
        j = 1; clear h;HYLIM = 0;
        for j = 1:length(UNIT)% plot psth for all unit
            %            clear h; % for  plot option 2
            %            HYLIM = 0; % for  plot option 2
            for icond = 1:nC
                if 1 % plot option 1: all conditions  on one subplot
                    h(j) = subplot(length(UNIT),1,j);
                    % % make plot pretty
                    sL{icond} = num2str(Val(icond,:));
                    stemp = sprintf('\\bf{Unit} %s\nSpk rate (Hz)',num2str(j));
                    if j == length(UNIT);                    xlabel('time (ms)'); end;
                else % plot option 2: each condition has its own subplot
                    sL = {};
                    h(icond) = subplot(nC,length(UNIT),j+(icond-1)*(length(UNIT)));
                    stemp = '';
                    if icond==1
                        title(['unit ' num2str(j) ' spks ' num2str(length(unit(j).sw))])
                        if j==length(UNIT);  stemp = 'Spk rate (Hz)'; else end
                        set(gca, 'XTick', [])
                    elseif icond==nC
                        xlabel('time (ms)');
                    else
                        set(gca, 'XTick', [])
                    end
                    stemp = sprintf('%s \n %1.2f',stemp, single(Val(icond,:)));
                end
                
                ylabel(stemp,'color',Unitcmap(j,:),'interpreter','Latex'); % color yaxis acording to units color
                
                if icond==1 & j ==1;    title(['PSTH ' sDesc]); end
                % % make plot pretty ends
                
                a = squeeze(PSTHcond(j,icond,:)); % mean
                xtime = [1:length(a)]*dt*1000;
                hold on
                try
                    b = filterdata(single(a),dt,5,0);%          b = smooth(a,.01/dt);
                catch ME % BA kluge
                    getReport(ME)
                    b = a;
                end
                plot(xtime,b/dt/PSTHnSw(j,icond),'Color',cmap(icond,:),'linewidth',2); plotset(1);
                axis tight;
                HYLIM = max(max(ylim),HYLIM);
            end
            if exist('h','var')
                linkaxes(h,'xy');axis tight;
                ylim([0 HYLIM]);
            end
            
        end
        % if bAnalyzeAllFilesTogether;
        sannotate = helpAnnotate(1,fid,FigureName,SParam);  %NESTED function (shares same workspace)
        
        % end
        if exist('sL','var'); if ~isempty(sL);            legend(sL);        end ; end
        
        %%% SUMMARY plots for different parameters
        if strfind(Vs,'Ori')  % this is a stupid,kloogy way of deal with different cases improve
            bdo = 1;
            if length(VarParam)>1
                bdo = 0;
                if isempty(strfind(VarParam(1).Name,'Contr')) & isempty(strfind(VarParam(2).Name,'Contr'))
                    bdo = 1;
                end
                if isempty(strfind(VarParam(1).Name,'Spa')) & isempty(strfind(VarParam(2).Name,'Spa'))
                    bdo = 1;
                end
            end
        end
        
        if isempty(analyzewindow) % use user specified window if not empty
            sampWIND = round([0 min(size(PSTHcond,3)-1,SParam.StimDuration/dt)]+1);
        else
            sampWIND = round([analyzewindow(1)/dt min(size(PSTHcond,3)-1,analyzewindow(2)/dt)]+1);
        end
        if bdo
            
            
            cmap = jetm(length(UNIT));
            fid = Nfig+20+F;figure(fid);set(fid,'Name', ['Orient Tuning ' sDesc]);clf;
            %             set(gcf,'Position',[2726+Nfig/2          -8         560
            %             420]);
            set(gcf,'Position',[ 2391+Nfig/2           33         560         420])
            clear Rpref Rortho;
            for j = 1:length(UNIT)
                if ~isempty(PSTHcond)
                    spikeInWindowPerSweep = (sum(PSTHcond(j,:,sampWIND(1):sampWIND(2)),3)./PSTHnSw(j,:))';
                    r = [spikeInWindowPerSweep; spikeInWindowPerSweep(1)]; theta = [Val; Val(1)];
                    try
                        subplot(3,1,1);
                        title(sDesc);
                        h = polar(deg2rad(theta),r/max(r)); hold all
                        set(h,'color',cmap(j,:),'linewidth',2)
                        subplot(3,1,2);
                        h2 = plot(theta(1:end-1),r(1:end-1),'-o'); hold all; plotset(1)
                        set(h2,'color',cmap(j,:),'linewidth',2)
                        temp = max(r(1:end-1)); % could be more than 1 max
                        indmax(1) = find(r(1:end-1)==temp(1),1,'first');% index of max response
                        indmax(2) = find(theta(1:end-1) == mod(theta(indmax(1))+180,360));% index of max response
                        if isempty(indmax(2));  indmax(2)  = NaN; end
                        
                        temp  = find(theta(1:end-1) == mod(theta(indmax(1))+90,360));
                        if isempty(temp); temp = NaN; end
                        indmaxortho(1) = temp ;
                        
                        temp  = find(theta(1:end-1) == mod(theta(indmax(1))-90,360));
                        if isempty(temp); temp = NaN; end
                        indmaxortho(2) = temp ;
                        
                        if ~isempty(indmaxortho) % can only calculate this if spiking at 90deg from max has been measured
                            Rpref(j) = mean(r(indmax(~isnan(indmax)))) ;% sum of respones to same orientation opposite direction
                            Rortho(j) = mean(r(indmaxortho(~isnan(indmaxortho))));
                        else
                            Rpref(j) = NaN;
                            Rortho(j)  = NaN;
                        end
                        
                        %** compute orienation selectivity using vector average
                        %(alternative method)
                        temp = r(1:end-1);
                        vavg = abs(temp.*exp(sqrt(-1)*mod(theta(1:end-1),180))/sum(temp));
                        % select perferred orientaiton and take mean (may be 2
                        % directions in PO
                        vavg = mean(vavg(~isnan(indmax)));
                        
                        if ~isnan(indmax(2)) ||~isnan(Rpref)
                            subplot(3,1,3);
                            % Note Cris does this by fitting 2 gaussians at theta
                            % and theata+pi first
                            % Orientation selectivity
                            if ~isnan(Rpref)
                                htemp(1) = plot(ones(size(Rpref(j))),(Rpref(j) - Rortho(j))./(Rpref(j) + Rortho(j)),'o');hold on;
                            end
                            % Dirn selectivity                                            % (Rpref-Roppsite)/(Rpref+Roppsite)
                            if ~isnan(indmax(2))
                                htemp(2)= plot(2*ones(size(r(indmax(1)))),(r(indmax(1)) - r(indmax(2)))./(r(indmax(1)) + r(indmax(2))),'o');
                            end
                            
                            htemp(3)= plot(3*ones(size(vavg)),vavg,'o'); % vector average
                            
                            set(htemp,'color',cmap(j,:));
                            set(gca,'XTick',[1 2 3],'XTickLabel',{'Osi','Dsi','Os_Vec'}); xlim([0.5 3.5]);
                            
                            
                            temp = ylim; ylim([0 max(1,temp(2))])
                            plotset(1)
                        end
                    catch % for debugging
                        icond
                    end
                    
                end
                subplot(3,1,2);
                ylabel('# spikes/sweep');
                set(gca,'XTick',[0:45:360])
                %             legend(sL);
                plotset(1);
                sannotate = helpAnnotate(1,fid,FigureName,SParam);
                
                subplot(3,1,3);
                ylabel('$$\frac{R_p- R_o}{R_p+ R_o}$$','Interpreter','latex');
                plotset(1)
            end
        elseif 1 %% size(Val,2)==1 % works for contrast, do nothing if analyzing all conditions
            %                         fid = Nfig+20+F;figure(fid);set(fid,'Name', [sDesc]);          clf;
            fid = 20+F;figure(fid);set(fid,'Name', [sDesc]);
            set(gcf,'Position',[1307         126         349         852]);
            %             all on one plot
            %             fid = 200+20+F;figure(fid); %set(fid,'Name', [sDesc]);
            
            %             set(gcf,'Position',[2726+Nfig/2          -8         560         420]);
            
            
            for i = 1:length(UNIT)
                if ~isempty(PSTHcond)
                    spikeRate = sum(PSTHcond(i,:,sampWIND(1):sampWIND(2)),3)./PSTHnSw(i,:)/(diff(sampWIND)*dt);
                    xval = Val;
                    bplotfit = 0;
                    if strfind(Vs,'Contrast') % fit
                        bplotfit = 1;
                        Starting=[max(y) 20 2 8]; % initial guess
                        % R = Rmax c^n/(c50^n + c^n) + B;
                        s = fitoptions('Method','NonlinearLeastSquares',...
                            'Lower',[0,0,0,0],...
                            'Upper',[1 1 1 1]*max(y)*5,...
                            'Startpoint',Starting);
                        fR = fittype('Rmax*(x^n/(c50^n + x^n)) + B','options',s);
                        [c2,gof2] = fit(x,y,fR)
                    end
                    
                    if AnalType == 2||AnalType == 1
                        cmap = jetm(NVAR); % color by Variable 1
                        
                        subplot(length(UNIT),1,i)
                        h2 = plot(xval,spikeRate,'-o');hold all;
                        set(h2,'color',cmap(kk,:),'linewidth',2);
                        %                     set(h2,{'DisplayName'},{num2str(VarParam(1).Values(kk
                        %                     ))})
                        set(h2,{'DisplayName'},{sDesc})
                        
                        if i ==length(UNIT); xlabel(VarParam(2).Name(1:3)); end
                        stemp = sprintf('Unit %s \n', num2str(i) );
                        ylabel([stemp 'Sp rate (Hz)']);
                        
                        
                    else
                        cmap = jetm(length(UNIT)); % color by unit
                        
                        h2 = plot(xval,spikeRate,'-o');hold on;
                        set(h2,'color',cmap(i,:),'linewidth',2)
                        set(h2,{'DisplayName'},{sDesc})
                        xlabel(Vs);
                        ylabel(['Sp rate (Hz)']);
                    end
                    if bplotfit
                        plot(c2,'r') ;
                        %                    XLIM = xlim; YLIM = ylim;
                        %                    stemp = sprintf('c_{50} = %1.1f',)
                        %                    text(max(XLIM)*.90 ,min(YLIM)*1.1)
                        %                    % ADD TITLE wit C50
                    end
                    hold on;
                end
                sannotate = helpAnnotate(1,fid,FigureName,SParam);  %NESTED function (shares same workspace)            title([num2str(WIND) ' sec']);
                %             YLIM = max(YLIM,ylim);
                %              ylim([0 YLIM])
                legend hide;legend show; legend boxoff
            end
        end
        %orientation tuning % for moving bar
        
    end

end
% ********************functions that don't share workspace
function [waves Nspks NBAD hsubplot] = spikequality(unitid,spkSortFN,hax)
% hax is optional, if given it should be a vector of 2 handles;
defpos(1,:) = [0 0.85 .1 .1];
defpos(2,:) = [0.015 0.70 .075 .1];

if nargin < 3
    hax(1) = axes('Position',defpos(1,:));
    hax(2) = axes('Position',defpos(2,:));
end

S_RIGSPECIFIC_SPIKESORT_CONFIG;
ABSREF = 3;
MAXSPIKESPLOT =500;

% ISI plot
smallWindow = 0.01; % sec

load(fullfile(SORTEDSPIKES_SAVEPATH,spkSortFN),'spikes','assign')
spkind = find(assign==unitid);  % plot only subset of spikes if there are too many

dt = 1/spikes.Fs;
Unitcmap = jetm(length(unique(assign(assign>0))));
tmp = max(floor(length(spkind)/MAXSPIKESPLOT)+1);
waves = spikes.waveforms(spkind(1:tmp:end),:);
orderstimes = sort(spikes.spiketimes(spkind));
isis =sort(diff(orderstimes));
NBAD =  sum(isis*1000<ABSREF);
Nspks = length(spkind);

t = dt*[1:size(waves,2)]*1e3; % ms
set(gcf,'CurrentAxes',hax(1));
% OVERLAYED PLOT
if ~isempty(waves)
    lh = mplot(t, waves, 'Color', Unitcmap(unitid,:)); axis tight; plotset(1);
    stemp = sprintf('Unit %d, Nspk: %d Nbad: %d',unitid,Nspks,NBAD); title(stemp);
    set(get(hax(1),'Title'),'color',Unitcmap(unitid,:));
    plotscalebar(4,{'mV','ms'});
end

% from chronux showclust
set(gcf,'CurrentAxes',hax(2));
isis = isis(isis <= smallWindow);
smalltimes = linspace(0,smallWindow,smallWindow*spikes.Fs);
if (~isempty(isis)), n = histc(isis,smalltimes);  else,  n = zeros(length(smalltimes));  end;
% lh = patch([0 ABSREF ABSREF 0]', [0 0 max(n) max(n)]', [0.1 0.1 0.1 0.1], [1 1 1], 'EdgeColor', [0.25 0.25 0.25],'FaceColor','none'); hold on;
lh = line([ABSREF ABSREF], [0 max(n)],'linewidth',2,'color','r'); hold on;
hp = plot(smalltimes*1000,n);    ylim = get(gca, 'YLim');
uistack(hp,'top'); % send patch the back
% set(hp,'FaceAlpha',0.5); % CAN'T do thos for some reason changes ALL
% other plots in the figure
set(gca,'Xlim',[0 smallWindow*1000]);
xlabel('ISI (ms)'); box off; plotset(1);
end
function maskout = maskspikes(iUnit,iCond,spikeData) % select spikes for analysis
% function maskout = maskspikes(iUnit,iCond,spikeData) % select spikes for analysis
% if iUnit or iCond are empty then ALL values are accepted
maskout = zeros(size(spikeData,1),1,'int16'); unitmask = ones(size(spikeData,1),1,'int16'); condmask = unitmask;
if ~isempty(iUnit)
    unitmask = ismember(spikeData(:,3),iUnit);
end
if ~isempty(iCond)
    condmask = ismember(spikeData(:,4),iCond);
end
maskout = unitmask&condmask;
% maskout = logical(maskout);
end



function [spikeData swcond Nspikes] = concatenateSpiketimesAcrossFiles(spikesortdata,indrawdatafilecond,rawdatainfo,includesweeps)
% function [spikeData swcond] = concatenateSpiketimesAcrossFiles(spikesortdata,indrawdatafilecond,rawdatainfo,includesweeps)
% NOTE: this function not only concatenates spikes from different files,
% but excludes spikes that occur outside of includedsweeps

DAQSAVEPATH = [];
S_RIGSPECIFIC_SPIKESORT_CONFIG;

% find spikes that are from sweeps that should be included.
for i = 1: length(indrawdatafilecond);
    indInclude(i) = {includeSweeps(spikesortdata.file(indrawdatafilecond(i)).sweep,includesweeps{i})};
end
Nspikes = 0;
% predefine array
for i = 1: length(indrawdatafilecond); Nspikes(i) = sum(indInclude{i}); end % length(spikesortdata.file(indrawdatafilecond(i)).spiketimes); end
spikeData = nan(sum(Nspikes),4);
temp1 = [0 cumsum(Nspikes)]+1;
swOffset = 0;
swcond = nan(1,cellfun(@sum,{rawdatainfo(indrawdatafilecond).maxtrigger})); %  predefine

for i = 1: length(indrawdatafilecond);
    tempind = find(indInclude{i});
    spikeData(temp1(i):temp1(i+1)-1,:) = [spikesortdata.file(indrawdatafilecond(i)).spiketimes(tempind)  ...
        (spikesortdata.file(indrawdatafilecond(i)).sweep(tempind) +swOffset) ...
        spikesortdata.file(indrawdatafilecond(i)).assign(tempind)  ...
        spikesortdata.file(indrawdatafilecond(i)).VstimCond(tempind) ];
    % sweeps are offset by number of sweeps in previous file other wise
    % spikes from the same sweep on different files will appear to happen
    % at the same time
    [junk temp] = getStimCond({fullfile(DAQSAVEPATH,rawdatainfo(indrawdatafilecond(i)).filename)});
    masknonnan = ~isnan(temp);
    stimmaxtrigger = sum(masknonnan);
    % check is same as DAQ maxtrigger
    if stimmaxtrigger~=rawdatainfo(indrawdatafilecond(i)).maxtrigger;
        % BA should fix DAQ and not allow this to happen
        printf('WATCH OUT !! StimCond file and DAQ file have different numbers of triggers');
        if (stimmaxtrigger-1)==rawdatainfo(indrawdatafilecond(i)).maxtrigger
            printf(' StimCondmaxtrigger = DAQmax trigger +1 so ASSUMED that LAST Stim Trigger does exist AND all others are OK');
            stimmaxtrigger = stimmaxtrigger-1;
        elseif rawdatainfo(indrawdatafilecond(i)).maxtrigger > (stimmaxtrigger)
            printf(' StimCondmaxtrigger < DAQmax trigger so ASSUMED unknown triggers set to NaN');
            stimmaxtrigger =rawdatainfo(indrawdatafilecond(i)).maxtrigger;
        else
            error('this difference was not fixed')
        end
        
    end
    swcond(swOffset+1:swOffset+rawdatainfo(indrawdatafilecond(i)).maxtrigger) = temp(1:stimmaxtrigger);
    
    swOffset = swOffset+rawdatainfo(indrawdatafilecond(i)).maxtrigger;
    
end
    function indInclude = includeSweeps(data,includesweeps)
        % find spikes that should be included
        if includesweeps
            indInclude  = ismember(data,includesweeps);
        elseif isempty(includesweeps)
            indInclude = ones(size(data));
        elseif includesweeps==0
            indInclude = zeros(size(data));
        end
        
    end
end
function [VarParam SParam] = formatVStimCellInfo(Vstimcellout,indSTFfilecond) % as struct for easy usage
rowN = indSTFfilecond; % remember must add one to STF fileind before it is passed in because the Vstimcell has titles in first row
% var params
if isnan(Vstimcellout{rowN,4}); NVar = 1; else NVar=2; end
for i = 1:NVar
    VarParam(i).Val= Vstimcellout{rowN,3+(i-1)*2};
    VarParam(i).Name = Vstimcellout{rowN,2+(i-1)*2};
end

%static parameters
SParam.StimName = Vstimcellout{rowN,6};
if strcmp( SParam.StimName , 'Drift Gratings') && Vstimcellout{rowN,21}
    SParam.StimName = 'DriftSq Gratings';
end
SParam.StimDuration = Vstimcellout{rowN,7};
SParam.Baseline = Vstimcellout{rowN,8};
SParam.tempFreq = Vstimcellout{rowN,12};

SParam.spfFreq = Vstimcellout{rowN,11};
SParam.bMask = Vstimcellout{rowN,29};

% String with SELECTED Extra bits of  information about stimulus
SParam.sDesc = '';
SParam.sDesc = sprintf('%s Mask %s\n', SParam.sDesc,num2str(Vstimcellout{rowN,9}));
SParam.sDesc = sprintf('%sSf %1.1f cpd\n', SParam.sDesc, Vstimcellout{rowN,11});
SParam.sDesc = sprintf('%sTf %1.1f Hz\n', SParam.sDesc, Vstimcellout{rowN,12});

temp = Vstimcellout{rowN,13}; % orienation
if max(size(unique(temp)))==1
    SParam.sDesc = sprintf('%s%1.0f Deg\n', SParam.sDesc, temp(1));
end

end

function [indexCol allcondname condVal] = getCond(cond,Vstimcellout)
% defines what cond.struct should contain
% and sets up condVal with either the values defined by the user in cond or
% default nan
condVal=struct([]);
indexCol = [29 9 14 11 12 21];% index of Col with condition in Vstimcellout
% NOTE: this should not overlap with indreq or indReqByFirst specified
% by the AnalType
allcondname = {'bMask','Mask','contrast','spfreq','tempFreq','squaregratings'};
for i = 1:size(allcondname,2); % which conditions are specified
    % check that column names are correct(if Vstimcellout changes format this will be screwed up)
    if ~isequal(allcondname{i},Vstimcellout{1,indexCol(i)});error('Vstimcellout Col names do NOT match conditions. Has Vstimcellout been changed?'); end
    condVal = checkcond(cond,lower(allcondname{i}),condVal); % get condVal
end

    function condVal = checkcond(cond,fieldn,condVal)
        % check cond struct, if a field exists get its value otherwise set its
        % value to nan
        if isempty(condVal)
            condVal = struct('tempfield',1); %(kluugy) add tempfield to get around fact that can't add a dynamic field to a struct that has not fields to begin with
        end
        if isfield(cond,fieldn)
            if isempty(cond.(fieldn)); % empty field is like not specifying fields
                condVal.(fieldn) = NaN;
            else
                condVal.(fieldn) = cond.(fieldn);
            end
        else
            condVal.(fieldn) = NaN;
        end
        if isfield(condVal,'tempfield')
            condVal = rmfield(condVal,'tempfield'); % remove temp field
        end
    end
end
function [indrawdatafilecond lstrawfilename ConditionDesc indSTFfilecond] = getfiles(cond,indreq,indReqByFirstFile,typereqcol,Vstimcellout,STF,rawdatainfo)
% function [fileindrawdata lstrawfilename ConditionDesc indSTFfilecond] = getfiles(cond,indreq,indReqByFirstFile,typereqcol,Vstimcellout,STF,rawdatainfo)
% is function gets the indexes (in rawdatainfo) of the files to analyze,
% this is used when extracting spikes from spikesortfile, to insure that
% only the spikes from the filenames that are correct for this analysis (AnalType" and the
% "cond" struct) are% extracted
%
% ConditionDesc is a cell array of the conditions that were met, (first row is column names) in all the
% rawdata selected for analysis. If there is more than one row in
% ConditionDesc rawdata that many ANY of the rows is included

[indexCol allcondname condVal] = getCond(cond,Vstimcellout);
filemask = zeros(1,size(Vstimcellout,1)-1);

ConditionDesc = cell(size(typereqcol,1)+1,size(Vstimcellout,2)); % predefine
ConditionDesc(1,:) = Vstimcellout(1,:);
% NOTE: may NOT be correctly dealing with conditions that take on 'all'

maskCondunspec = logical(cellfun(@any,cellfun(@isnan,struct2cell(condVal),'UniformOutput',false))); %     find  unspecified conditions
% find conditions NOT specifed as 'all' any value will be accepted for this
% condition % BA don't understand what I was doing here FIGURE OUT!
indNOTall = find(cellfun(@isempty,struct2cell(structfun(@(x) strfind('all',x),condVal,'UniformOutput',false))));
% remove  conditions specifeid as 'all'
tempcell = struct2cell(condVal)'; % cell of spec. cond
tempcell = tempcell(indNOTall);            indexCol = indexCol(indNOTall);            maskCondunspec = maskCondunspec(indNOTall);
% set nans to zero because nan~=nan
tempcell(cellfun(@any,cellfun(@isnan,tempcell,'UniformOutput',false))) = {0};
tempVstimcellout = Vstimcellout; tempVstimcellout(cellfun(@any,cellfun(@isnan,tempVstimcellout,'UniformOutput',false))) = {0};
typereqcol(cellfun(@any,cellfun(@isnan,typereqcol,'UniformOutput',false))) = {0};
indTemplateFile = 0; % index of first "Templatefile" that matches conditions and reqconditions first file
for i = 2: size(Vstimcellout,1) % find Templatefile
    REQ = 0;
    for j = 1 :size(typereqcol,1)
        REQ = REQ || isequal(typereqcol(j,:),tempVstimcellout(i,indreq(1,:)));
        ConditionDesc(j+1,indreq(1,:)) = typereqcol(j,:);
    end
    %     if size(typereqcol,1)== 2 % allowed to satisfy 1 of 2 conditions if typereqcol and indreq have 2 rows
    %         REQ = isequal(typereqcol(1,:),tempVstimcellout(i,indreq(1,:))) ...
    %             || isequal(typereqcol(2,:),tempVstimcellout(i,indreq(2,:)));
    %
    %         ConditionDesc(2,indreq(1,:)) = typereqcol(1,:);
    %         ConditionDesc(3,indreq(1,:)) = typereqcol(2,:);
    %
    %     else
    %         REQ = isequal(typereqcol(1,:),tempVstimcellout(i,indreq(1,:)));
    %         ConditionDesc(2,indreq(1,:)) = typereqcol(1,:);
    %
    %     end
    if ~indTemplateFile %  has Templatefile  been found?
        if REQ ... % required conditions (see above)
                &&     isequal(tempcell(~maskCondunspec),tempVstimcellout(i,indexCol(~maskCondunspec)))%  match specified conditions as well
            indTemplateFile = i; % this is templatefile
            filemask(i-1) = 1;
        end
    elseif REQ ... % match required, specified, and unspecified (to Template)
            &&     isequal(tempcell(~maskCondunspec),tempVstimcellout(i,indexCol(~maskCondunspec)))...
            &&     isequal(tempVstimcellout(indTemplateFile,[indReqByFirstFile indexCol(maskCondunspec)]),tempVstimcellout(i,[indReqByFirstFile indexCol(maskCondunspec)])) % for unspecified condition files must match condition of TemplateFile
        filemask(i-1) = 1;
    end
end
if indTemplateFile % if a file is found
    for j=2:size(ConditionDesc,1)
        ConditionDesc(j,[indReqByFirstFile indexCol(maskCondunspec)]) = tempVstimcellout(indTemplateFile,[indReqByFirstFile indexCol(maskCondunspec)]);
        ConditionDesc(j,indexCol(~maskCondunspec)) = tempcell(~maskCondunspec);
    end
end
% now find the index of files to analyze in 'spikesortdata' by looking at  'rawdatainfo'
clear rawfilenames;
for i =1:size(rawdatainfo,2)
    tempind  = strfind(rawdatainfo(i).filename,'.');
    if isempty(tempind); tempind = length(rawdatainfo(i).filename)+1; end
    rawfilenames{i} = rawdatainfo(i).filename(1:tempind-1);
end
indSTFfilecond = find(filemask); % index of files in STF that meet conditions
indrawdatafilecond = nan(length(indSTFfilecond),1); %index of STFfiles in spikesort file (and rawdatainfo)
%  STF doesn't have to specified in the same order as the rawfilenames were spikesorted
for i =1:length(indrawdatafilecond) % for all STF.filenames that met conditions
    indrawdatafilecond(i) = find(cellfun(@(x) strcmp(STF(indSTFfilecond(i)).filename,x),rawfilenames));
    % note Nans here mean raw daq files hasn't been analyzed
end
lstrawfilename = rawfilenames(indrawdatafilecond);

% check that file have same length and sampling rate
for i =2:length(indrawdatafilecond)
    if rawdatainfo(indrawdatafilecond(1)).dt ~= rawdatainfo(indrawdatafilecond(i)).dt ||...
            rawdatainfo(indrawdatafilecond(1)).maxsample ~= rawdatainfo(indrawdatafilecond(i)).maxsample
        error('Specified Filenames cannot be analyzed together they have different length or sampling rate');
    end
end

end


function [tuning indmax] = calcDirnOriTuning(r,theta)
%     function [tuning indmax] = calcDirnOriTuning(r,theta)
% outputs struct with tuning parameters
% and ind of preferred direction in r

% calcualted Orientation/ Dirn SelectivityIndex (see Niell Stryker)
%         (Rpref(j) - Rortho(j))./(Rpref(j) + Rortho(j))
%** compute orienation selectivity using vector average

temp = max(r(1:end-1)); % could be more than 1 max
if ~isempty(temp)||~all(isnan(temp))
    % orientation (get both directns)
    indmax(1) = find(r(1:end-1)==temp(1),1,'first');% index of max response
    temp2 =  find(theta(1:end-1) == mod(theta(indmax(1))+180,360)); % index of 180 max response
    if isempty(temp2);  indmax(2)  = NaN; else indmax(2) = temp2; end % because can't set indmax(2) = []; for some reason
    %%get orthogonal directions
    temp  = find(theta(1:end-1) == mod(theta(indmax(1))+90,360));
    if isempty(temp); temp = NaN; end
    indmaxortho(1) = temp ;
    temp  = find(theta(1:end-1) == mod(theta(indmax(1))-90,360));
    if isempty(temp); temp = NaN; end
    indmaxortho(2) = temp ;
    
    if ~isempty(indmaxortho) % can only calculate this if spiking at 90deg from max has been measured
        Rpref = mean(r(indmax(~isnan(indmax)))) ;% sum of respones to same orientation opposite direction
        Rortho = mean(r(indmaxortho(~isnan(indmaxortho))));
    else
        Rpref = NaN;
        Rortho  = NaN;
    end
    
    tuning.OSI = (Rpref - Rortho)./(Rpref + Rortho);
    if ~isnan(indmax(2))    tuning.DSI = (r(indmax(1)) - r(indmax(2)))./(r(indmax(1)) + r(indmax(2)));
    else tuning.DSI = nan; end
    
    tuning.prefD = mod(theta(indmax(1))+90,360);
    
    %(alternative vector average method)
    temp = r(1:end-1);
    VAVG = abs(temp.*exp(sqrt(-1)*mod(theta(1:end-1),180))/sum(temp));
    % select perferred orientaiton and take mean (may be 2
    % directions in PO
    tuning.VAVG = mean(VAVG(~isnan(indmax)));
else
    tuning.OSI = NaN; tuning.DSI = NaN; tuning.VAVG = NaN;
end

end
function F1F0 = getF1F0(y,dt,stimfreq)
res = 0.5; % Hz
% set freq resolution to res Hz unless trace is too short otherwise
% max allowed
nfft = 2^nextpow2(min((1/res)/dt,length(y)));
Y = fft(single(y),nfft);
Pyy = Y.* conj(Y) / nfft;
Pyy = Pyy(1:nfft/2+1);
% Graph the first 514 points (the other 512 points are redundant) on a meaningful frequency axis:
f = 1/dt*(0:nfft/2)/nfft;
temp = abs(f-stimfreq);
ind = find(temp==min(temp),1,'first');
F1=Pyy(ind);
F0 = Pyy(1);

F1F0 = F1/F0;
if 0 % debug
    figure(1000);clf;
    plot(f,Pyy,'-r')
    xlim([0 10])
    title('Frequency content of y')
    xlabel('frequency (Hz)')
    %
    %  % alternative method
    % param.tapers = 1;
    % param.Fs=1/dt;
    % param.fpass = [0 20]
    % param.fpass
    % y = squeeze(PSTHcond(i,indmax(1),sampWIND(1):sampWIND(2)));
    % plotpsCHX(y,dt,2000,param.fpass,param,[],1)
end
end
function [c2 gof2]=fitcontrast(x,y)
% fit hyperbolic ratio
Starting=[max(y) 0.2 2 min(y)]; % initial guess

i = 1; gof2.rsquare = 0;
while (i<5 && gof2.rsquare<0.6)  % until a fit is found
    % R = Rmax c^n/(c50^n + c^n) + B;
    s = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[0,0,0,0],...
        'Upper',[1 1 1 1]*max(y)*1.5,...
        'Startpoint',Starting);
    fR = fittype('Rmax*(x^n/(c50^n + x^n)) + B','options',s);
    [c2,gof2] = fit(x,y,fR);
    i =  i+1;
    Starting(2:end)=Starting(2:end)*rand(1); % switch initial guess
end
if gof2.rsquare<0.8
    warning('**************** Bad CRF fit')
end
end
function setColorbyColorOrder(hplot,hline,nplots)
% function 1) gets the colororder index (PlotColorIndex) for the next line from hplot(1)
%     (which would autoincrement
% everytime you plot if 'hold all' is set).
% 2) sets the hline handles to this color
% 3) increments the colororder index in all hplots
% INPUT
%    hplot(handles for axes/subplots), should be a vector
%     hline (handles for lines) can be a cell row vector all handles will be concatentated and
%     manimupated together
%    nplots (optional) if more than 1 plot is being plotted on the axis
%    and you want N of these plots to augment the colororder
%    by just 1, specify the TOTAL number of plots here (otherwise will
%    augment  by number of plots)
% note - must use 'hold all'
%

if nargin<3; nplots = 1; end;
cind = getappdata(hplot(1),'PlotColorIndex'); % assumes all the plots have same state as far as the number lines plotted on them
holdmode = getappdata(hplot(1),'PlotHoldStyle'); % assumes all the plots have same state as far as the number lines plotted on them
Co = get(hplot(1),'colorOrder');
cind = cind-(nplots-1);
if isempty(cind) || cind<=1 % first plot
    cind = 1;
    helpsetcolor(hline,Co(cind,:));
    VECTsetappdata(hplot,'PlotColorIndex',cind+1);
    
elseif cind >=2;
    if holdmode % % hold all
        helpsetcolor(hline,Co(cind-1,:));
    else
        error('must use hold all')
    end
    VECTsetappdata(hplot,'PlotColorIndex',cind);
    
end


end

function VECTsetappdata(H, NAME, VALUE)
%     function VECTsetappdata(H, NAME, VALUE)
% like setappdata fun but H can be a vector not just a scalar.
% all H are set to she same VALUE
for i = 1:length(H)
    setappdata(H(i),NAME,VALUE);
end

end

function helpsetcolor(ph,C)

h = [];
for i = 1:size(ph,2)
    if ~isempty(ph{i})
        h = [h; ph{i}];
    end
end
if ~isempty(h)
    set(h,'Color',C);
end
end



function switched = helpswcond(data,nV2,nV1)
% for switching VarCond1 and VarCond2
V1 = ceil(data/nV2);
V2 = data - (V1-1)*nV2;
% switched = (V2-1)*nV2+V1;
switched =  (V2-1)*nV1+V1;
% SHOULD BE
%  (V2-1)*nV1+V1;
% otherwise too big (V2-1) *nV2
% BUT why doesn't it work for other files?
% [V1 V2 data switched]
end

function savetuningProperties(newdata,newdataCond)
%           newdata and newdataCond are structs.
%           newdata must contain the following fields.
%           They are used to make each unit unique, and check whether the same unit exists already in the struct
%           (see below for more details)
%           newdata.spkSortFN % filename of spikesortfile of unit
%                  .Unit % unit number in spikesortfile assign
%                  .rawdata_in_analysis % rawdata files includes in analysis
%
% - A unit(unique by Unit and sortfile entry) can appear multiple times in
% this table because if you ANALYZE the same unit with data from different
% files. While this is useful for comparing different conditions it could
% be confusing.
%  - The PSC Vstim conditions (and others?) are defined by the
% rawfilename and not in this table (but can be extracted by using the
% rawfilename to look up where they are stored (e.g using  updateVStimCondTable)
% - UnitISOLATION information and spikewaveform, spiketime, information is also not saved here but
% can be retrieaved from the spksortfile.
% - Some of the analyzed data is sorted here e.g. PSTH and average tunig
% curve. These data is somewhat redundent, but I thought it would be
% convient not to have to regenerate.
% - A path to the figures created when the data was analyzed is included
%    all headers of this kind start with PATH_

S_RIGSPECIFIC_SPIKESORT_CONFIG;
filename = 'UnitTuningStruct';
try
    load(fullfile(TABLES_SAVEPATH,filename));
    
    % check if potentially overlaping entries already exist:
    % if same , exptname,  sortfile,  unit number, and the
    % ANY overlap with same rawdata files exist, replace the entry
    % NOTE:
    %     1)if analysis were setup to analyze part of one file at one
    % time and part in another then this approach would not
    % allow
    % the data to be both saved
    %     2) doesn't protect you from including data twice resorting
    %     (because a different unit name could result or a different
    %     sortfile name)
    
    % delete all rows that have same spkSortFn, unit
    % number, and same rawdata files (can have same unit but
    % different conditions in different rawdata files).
    % NOTE : a unit must have exactly the same rawdata
    % files to be overwritten. unit tuning with only a
    % partial overlap or not removed
    
    % mask of entries with same spkSortFN
    %                     indssfn = ~cellfun(@isempty,cellfun(@(x) findstr(FilenameSpkSort,x),UtuningTable(:,5),'UniformOutput',false));
    indssfn = ~cellfun(@isempty,cellfun(@(x) findstr(newdata.spkSortFN,x),{unitsummary.spkSortFN},'UniformOutput',false));
    % mask of entries with same unit number
    %                     indunit = ~cellfun(@isempty,cellfun(@(x) findstr(unt,x),UtuningTable(:,4),'UniformOutput',false));
    indunit = ~cellfun(@isempty,cellfun(@(x) findstr(newdata.Unit,x),{unitsummary.Unit},'UniformOutput',false));
    % should not be possilbe to have the same spikesortfile and
    % different exptname (so this is not checked)
    indtemp = find(indssfn & indunit); inddel = zeros(size(indunit));
    
    for ii = 1:length(indtemp) % check if rawdatafiles analyzed are EXACTLY the same as existing unit and file ( if so delete
        if isequal(unitsummary(indtemp(ii)).rawdata_in_analysis,newdata.rawdata_in_analysis)
            %                         mark for deletion (don't delete here because
            %                         changes indices of cell when a row is deleted.
            inddel(indtemp(ii)) = 1;
        end
    end
    % remove inddel entries
    unitsummary = orderfields(unitsummary(~inddel)); % note this seems like a slow way to delete because the entire struct must be recopied (should be a better way)
    unitsummaryCond = orderfields(unitsummaryCond(~inddel));
    ind = length(unitsummary)+1;
    
catch % create file if non exists
    ind =1;
end
%  add new entry
ftemp = fieldnames(newdata) ; for i = 1:length(ftemp); unitsummary(ind).(ftemp{i})=newdata.(ftemp{i});  end % must add field names one by one, cause cell2struct doesn't work if struct doesn't exactly match existing strtuct
ftemp = fieldnames(newdataCond) ; for i = 1:length(ftemp); unitsummaryCond(ind).(ftemp{i})=newdataCond.(ftemp{i}); end %
% unitsummary(ind) = cell2struct(struct2cell(newdata),fieldnames(newdata));

save(fullfile(TABLES_SAVEPATH,filename),'unitsummary','unitsummaryCond');
end