function summarizeTuningofUnits(Units,FilenameSpkSort,lstrawfilename,fid)
% function summarizeTuningofUnits(Units,FilenameSpkSort,lstrawfilename,fid)
%
% USAGE: summarizes Orientation tuning of units in the same spikesortfile,
% with rawdatafiles = lstrawfilename
%
% lstrawfilename is a cell array of rawfilenames. unit tuning data must be collected
% from filenames that are ALL within lstrawfilename. this allows the function to be used
% to summarize only cells that meet specific conditions (use function
% getfiles to determin the lstrawefilename for a certain set of condtions)
% load all or subset of units associated with FilnameSpkSort
% plot summary of their tuning properties
% % SUMMARY OF UNITS
if ~iscell(lstrawfilename); lstrawfilename= {lstrawfilename}; end

s_RIGSPECIFIC_SPIKESORT_CONFIG;

nUnitsSumarized = length(Units);

% plot waveforms
load(fullfile(SORTEDSPIKES_SAVEPATH,FilenameSpkSort),'spikes','assign')
filename = 'UnitTuningStruct';
load(fullfile(TABLES_SAVEPATH,filename));       ifidlist = []; % list of figure handles
for in = 1:length(Units)
    % mask of entries with same spksortfn
    indssfn = ~cellfun(@isempty,cellfun(@(x) findstr(FilenameSpkSort,x),{unitsummary.spkSortFN},'UniformOutput',false));
    % mask of entries with same unit number
    indunit = ~cellfun(@isempty,cellfun(@(x) findstr(Units(in),x),{unitsummary.Unit},'UniformOutput',false));
    % load Vstim params
    uind = find(indssfn& indunit);
    
    % get list of rawdatafiles that meet conditions
    tempind = zeros(size(uind));
    for k =1 : length(uind)
        if all(ismember(unitsummary(uind(k)).rawdata_in_analysis,lstrawfilename)) ; % rawdatafiles included in Unit analysis must ALL be within lstrawfilenames
            tempind(k) = 1;
        end
    end
    uind = uind(logical(tempind));
    
    
    annot = cell(length(uind),1);
    for k = 1:length(uind) %more than one condition for a unit can result in more than one index
        % fix nRow and nCol
        sexptname = unitsummary(uind(k)).ExptName;
        nunit =unitsummary(uind(k)).Unit;
        rateMax = [unitsummary(uind(k)).max_rate_hz unitsummary(uind(k)).BL_max_rate];
        rateMed =[unitsummary(uind(k)).med_rate unitsummary(uind(k)).BL_med_rate];
        prefDirn =  unitsummary(uind(k)).pref_Dirn ;
        spfFreq = unitsummary(uind(k)).spf_Freq;
        tuntemp.OSI =    unitsummary(uind(k)).OSI ;
        tuntemp.DSI = unitsummary(uind(k)).DSI;
        tuntemp.VAVG =  unitsummary(uind(k)).Vavg;
        F1F0 = unitsummary(uind(k)).F1F0 ;
        PSTHcond = unitsummary(uind(k)).PSTH ;
        PSTHnSw =  unitsummary(uind(k)).PSTHnSw;
        maxsample = unitsummaryCond(uind(k)).maxsample;
        maxtrigger = unitsummaryCond(uind(k)).maxtrigger;
        dt = unitsummaryCond(uind(k)).dt ;
        CondVal =  unitsummaryCond(uind(k)).CondVal;
        sampWIND  = unitsummaryCond(uind(k)).sampWIND;
        baselinesampWIND = unitsummaryCond(uind(k)).baselinesampWIND;
        totalUNITs= unitsummaryCond(uind(k)).totalUNITs;
        tempann = unitsummaryCond(uind(k)).sannotate ;
        % if different Vstim details for different units then cover up
        % annotation (should fix plotAnn so can replace and not cover)
        ConditionDesc =  unitsummaryCond(uind(k)).ConditionDesc;
         % WHY ARE bMASK and not MASK plotted together?
        
            % check is just orientation is variing 
        if isequal('Orientation',ConditionDesc{2,2}) && ~isstr(ConditionDesc{2,4})
            % **make sure that same conditions are plotted on same plot
            if ~isempty(ifidlist) % look through figures and find the one with the same conditions (e.g. mask, orientation etc specified when unit was analyzed)
                j = 1;
                while(~isequal(ifidCond{j},ConditionDesc))
                    j = j+1;
                    if j>length(ifidlist) % no similar figure exists make a new figure
                        ifidCond{end+1} = ConditionDesc;
                        ifidlist(end+1) =  ifidlist(end)+1; 
                    end
                end
                 ifid = ifidlist(j);
            else ifid = fid + 1; ifidlist(1) = ifid; ifidCond{1} =  ConditionDesc; end
                           

%             ifid
            figure(ifid); orient landscape;
            
            
            Unitcmap =  jetm(length(totalUNITs)); % color units
            
            
%             if exist('annot','var') && ~isempty(annot{k}); if ~isequal(regexprep(annot{k},'\s',''),regexprep(tempann{k},'\s',''));
%                     annot{k} = sprintf('          \n          \n          \n          \n'); end
%             else annot{k} = tempann{k}; end
            plotAnn( regexprep(tempann,'\s',' '),ifid,1)
            
            try
                bskipPSTH = 0;
                if size(PSTHcond,1)==nCond(ifid); %define figure parameters
                else display('Skipped PSTH'); bskipPSTH = 1; end% bskipPSTH (must skip psth cause different number of conditions for this unit
            catch % case where there is nCond is not defined yet
                nCond(ifid) = size(PSTHcond,1); % assum n conditions is same for all files
                nCol(ifid)= max([4,nUnitsSumarized,nCond(ifid)]); % at least 4 colns other wise the bigger of the number of units and num of conditions
            end
            nRow = 4;
            
            spkind = find(assign==nunit); MAXSPIKESPLOT =500; % plot only subset of spikes if there are too many
            tmp = max(floor(length(spkind)/MAXSPIKESPLOT)+1);
            t = dt*[1:size(spikes.waveforms,2)]*1e3; % ms
            
            hs(in) = subplot(nRow,nCol(ifid),[1+(in-1)*floor(nCol(ifid)/nUnitsSumarized):(in)*floor(nCol(ifid)/nUnitsSumarized)]);
            waves = spikes.waveforms(spkind(1:tmp:end),:);
            % OVERLAYED PLOT
            if ~isempty(waves)
                lh = mplot(t, waves, 'Color', Unitcmap(Units(in),:)); axis tight; plotset(1);
                stemp = sprintf('Unit %d, Nspk: %d',Units(in),length(spkind)); title(stemp);
                set(get(hs(in),'Title'),'color',Unitcmap(Units(in),:));
                if in==1; plotscalebar(4,{'mV','ms'}); else axis off; end
            end
            
            % legend
            sl{in} = [regexprep(sexptname,'\s',' ') 'U: ' num2str(Units(in))];
            
            xtime = [1:size(PSTHcond,2)]*dt;
            %%% PSTH
            if ~ bskipPSTH % skip when nCond doesn't match nCond of previously plotted unit
                for j = 1: nCond(ifid) %NOTE This PSTH plot only works if all units have same number of conditions (it is probalby technically not the best tocompare other metrics like OSI and DSI when units were probed with different number of angles..
                    
                    clear smoothPSTHcond;
                    if sum(PSTHcond(j,:))>0 % save time
                        smoothPSTHcond(j,:) = filterdata(single(squeeze(PSTHcond(j,:))),dt,5,0);%          b = smooth(a,.01/dt);
                    else                             smoothPSTHcond(j,:) = nan(1,size(PSTHcond(j,:),2)); end
                    
                    %%% plot
                    
                    h1(k,j) = subplot(nRow,nCol(ifid),[(nCol(ifid)+(j-1)*floor(nCol(ifid)/nCond(ifid))+1):(nCol(ifid)+(j)*floor(nCol(ifid)/nCond(ifid)))]);
                    plot(xtime,squeeze(smoothPSTHcond(j,:))/dt/PSTHnSw(j),'linewidth',2,'color',Unitcmap(Units(in),:)); hold all; plotset(1);
                    axis tight;
                    % baseline period baselinesampWIND
                    line(baselinesampWIND(1)*dt.*[1 1],[0 5],'color',[0 0 0])
                    title(num2str(CondVal(j)));
                end
                if ~isempty(h1); %                     format axis
                    stemp = sprintf('spk rate (Hz) **CHECK');
                    set(get(h1(k,1),'YLabel'),'String',stemp);
                    set(get(h1(k,1),'XLabel'),'String','sec');
                    temphh = h1(k,:); temphh = temphh(temphh>0);
                    set(temphh(max(1,end-3)),'XTickLabel',[],'YTickLabel',[]);
                    temp = ([min([cellfun(@min,get(temphh,'YLIM'));ylim']) max([cellfun(@max,get(temphh,'YLIM'));ylim'])]);
                    linkprop(temphh,{'box','Color','XLIM'});
                    xlim([ 0 maxsample*dt]);set(temphh,'YLIM',temp);
                end
            end
            % Summarize across tuning properties
            
            % plotting tunign properties
            
            spikeRate = sum(PSTHcond(:,sampWIND(1):sampWIND(2)),2)'./PSTHnSw/(diff(sampWIND)*dt);
            r = [spikeRate'; spikeRate(1)]; r(isnan(r)) = 0; %replace nan with zero WATCH out this could be bad..
            theta = [CondVal'; CondVal(1)];
            if ~isnan(baselinesampWIND) % baseline calculated for all conditions
                baselinespikeRate = sum(PSTHcond(:,baselinesampWIND(1):baselinesampWIND(2)),2)'./PSTHnSw/(diff(baselinesampWIND)*dt);
                rBL = [baselinespikeRate'; baselinespikeRate(1)]; rBL(isnan(rBL)) = 0;
            else                       rBL = nan(size(theta));                    end
            
            subplot(nRow,nCol(ifid),nCol(ifid)*2.*[1 1]+[1 nCol(ifid)/2]); % subplot takes half the row % polar plot
            h2(1,in,k) = plot(theta(1:end-1),r(1:end-1),'-o'); hold all; plotset(1)
            set(h2(1,in,k),'linewidth',2,'color',Unitcmap(Units(in),:));
            h2(2,in,k) = plot(theta(1:end-1),rBL(1:end-1),'-');  plotset(1) % baseline
            set(h2(2,in,k),'color',Unitcmap(Units(in),:));
            ylabel('rate (Hz)');
            set(gca,'XTick',[0:45:360])
            
            hpolar = subplot(nRow,nCol(ifid),nCol(ifid)*2.*[1 1]+[nCol(ifid)/2+1 nCol(ifid)]); % subplot takes the next half
            h3(1,in,k) = polar(deg2rad(theta),r/max(r)); hold all
            set( h3(:,in,k),'linewidth',2,'color',Unitcmap(Units(in),:))
            
            subplot(nRow,nCol(ifid),nCol(ifid)*3+[1:floor(nCol(ifid)/2)]);
            h4(1,in,k)= plot(ones(size(tuntemp)),tuntemp.OSI,'o');hold all;
            % Dirn selectivity                                            % (Rpref-Roppsite)/(Rpref+Roppsite)
            h4(2,in,k)= plot(2*ones(size(tuntemp)),tuntemp.DSI,'o');
            h4(3,in,k)= plot(3*ones(size(tuntemp)),tuntemp.VAVG,'o'); % vector average
            h4(4,in,k)= plot(4*ones(size(F1F0)),F1F0,'o'); % vector average
            set(gca,'XTick',[1:4],'XTickLabel',{'Osi','Dsi','Os_Vec','F1F0'}); xlim([0.5 4+.5]);
            temp = ylim; ylim([0 1])
            ylabel('$$\frac{R_p- R_o}{R_p+ R_o}$$','Interpreter','latex');   plotset(1)
            set(h4(:,in,k),'color',Unitcmap(Units(in),:));
            
            subplot(nRow,nCol(ifid),nCol(ifid)*3+[floor(nCol(ifid)/2)+1:nCol(ifid)]);
            indexX = 1;
            h5(1,in,k)= plot(1,rateMax(1),'o'); hold all; % stimulus
            h5(2,in,k)= plot(2,rateMax(2),'x');% baseline
            h5(3,in,k)= plot(3,rateMed(1),'o'); % stimulus
            h5(4,in,k)= plot(4,rateMed(2),'x');% baseline
            indexX = 4;
            
            set(gca,'XTick',[1:indexX],'XTickLabel',{'Max','BLMax','Med' ,'BLMed' }); axis tight; xlim([0.5 indexX+.5]);
            ylabel('rate(Hz)');   plotset(1);
            set(h5(:,in,k),'color',Unitcmap(Units(in),:));
            
            
            
        end
    end
    if length(uind)
        hs = hs(hs>0);
        if length(hs)>1
            temp = ([min(cellfun(@min,get(hs,'YLIM'))) max(cellfun(@max,get(hs,'YLIM')))]);
            set(hs,'YLIM',temp);
        end
        linkprop(hs,{'box','Color','YLIM','XLIM'});
    end
    
end
subplot(hpolar)
sl(cellfun(@isempty,sl)) = ''; % deal with empty strings

legend(sl,'box','off','color','none','position',[0.5 0.437 0.143 0.072])

end

