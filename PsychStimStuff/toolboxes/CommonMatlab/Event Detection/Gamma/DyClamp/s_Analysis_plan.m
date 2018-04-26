%% Algorithm for finding average membrane potential triggered on 
% LFP (or whole cell current).
% CASE 1: conductance was injected
% CASE 2: No conductance injected.

% algorigh steps
% load daata
output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_10_25_0008_1-2.abf',-1,1);
%% ORDER of channels in 
%% intracellular
%% Injected Currents
%% Monitor of GAmma (Intra or LFP)
indexoffset = 1; b_intra= 0;
%%Monitor
intradata = double(output.data(:,indexoffset+1:output.Nchan:end))*output.gain(1) + output.offset(1);
% remove mean of each sweep
intradata = intradata - repmat(mean(intradata),length(intradata),1);
curdata = output.data(:,indexoffset+2:output.Nchan:end);
mondata = double(output.data(:,indexoffset+5:output.Nchan:end))*output.gain(3);
mondata = (mondata-mean(mean(mondata)));

%% low pass filter MONITOR data
LPFilter = 60 ;%% low pass filter of LFP
ISIrange = [1/15 1/60]*1000; % ms
[B,A] = butter(2,2*LPFilter*output.dt,'low');
prodata =   filtfilt(B,A,mondata);

figure(1) %% PLOT to check threshold
clf
plot(prodata(1:20000,2))
title('Click to SET threshold')
[x y] = ginput(1)
thres = y;

findOSC_detectOSC002 %% much faster than 001

% FIND Max of Membrane Potential in window around LFP alignment



%%% CREATE REPORT:
%% look at distribution of "rise time of LFP"
% indEvent =Strough; % temporary should be
indEvent = double(eventsw).*length(prodata)+timeinsw;

%% TIME of injected current.
% remove mode of trace (assume mostly at baseline)
[a x] = hist(double(curdata(:,1)),500);
temp = max(a); temp = x(find(a == temp(1)));
curdata = curdata - int16(temp)*ones(size(curdata),'int16'); % subtract baseline
threscur = -1000; %% USER DEFINED THRESHOLD



cinj_ind = find(diff(curdata < threscur) ==1)+1; %% Conductance injection TIME in experiment of 
intv = 200; %(ms) USER DEFINE interval that conductances should be injected at
syndur = 32;  %(ms) USER DEFINE duration of synaptic conductance
if max(diff(cinj_ind)< intv/1000/output.dt/2)%% there can be noise check that in correct crossings are not detected
    warning ('Current Injection times (cinj_ind) may be wrong. Too short an inteval detected');
end


%% find "LARGE" current injection.
 %% ASSUMES NEGATIVE CURRENT
cinj_indL = find(diff(curdata < -1*.8*max(abs(unique(curdata))))==1)+1;
temp = findSimEvent(cinj_ind,cinj_indL,round(32/1000/output.dt)); %% note this only works if intv betweenconductance and the duration of conductance don't overlap
 cinj_indL = cinj_ind(temp(:,1));
 temp2 = cinj_ind; temp2(temp(:,1)) = 0;
 cinj_indS = nonzeros(temp2);

sprintf('REPORT// Total ConInj: %d\tLarge: %d\tSmall: %d', length(cinj_ind),length(cinj_indL),length(cinj_indS) )
if (length(cinj_indS) + length(cinj_indL) ~=  length(cinj_ind))
    error('cinj_small and cinj_lage don''t add up to Total cinj_ind')
end


%% find APs
thresAP = -10; %(mV) %% USER DEFINE Action potential threshold
indAP = find(diff(intradata(:,1)>thresAP)==1);
if ~isempty(indAP)
    wind = 50 %(ms) ind_sid within +/- window of cinj_inL are classifed as occuring with cinj_inL
    temp1 = findSimEvent(indAP,indEvent,round(wind/1000/output.dt))
    evSup = indEvent(temp1(:,2)) % superthres events
    temp =  indEvent;%
    temp(temp1(:,2)) = 0;
    evSub = nonzeros(temp);
    sprintf('REPORT// GammaCyc: %d\tSup: %d\Sub: %d', length(indEvent),length(evSup),length(evSub) )
else
    evSub = indEvent;
    evSup = [];
end

%% CLASSIFICATION of GAMMA events:
%     cinj subCTR  cinjsub cinjsubS cinjsubL  supCTR  cinjsup cinjsupS cinjsupL 
% COL   1    2        3         4        5         6     7      8      9
gammaClass = zeros(length(indEvent),9,'int16');
cinj = 1; csubctr = 2; csub = 3;csubS = 4;csubL = 5;
          csupctr = 6; csup = 7;csupS = 8;csupL = 9;
wind = 50;
%% SUB threshold
% WOI around gamma oscillation e.g 50ms. (SAME AS ABOVE
temp = findSimEvent(cinj_ind,indEvent,round(wind/1000/output.dt));
gammaClass(temp(:,2),cinj)= 1; %% index of events with INJ
temp = findSimEvent(cinj_ind,evSub,round(wind/1000/output.dt));
gammaClass(temp(:,2),csub)= 1;  %% Sub
temp = findSimEvent(cinj_indS,evSub,round(wind/1000/output.dt));
gammaClass(temp(:,2),csubS)= 1;  %% SubSmall
temp = findSimEvent(cinj_indL,evSub,round(wind/1000/output.dt));
gammaClass(temp(:,2),csubL) = 1; %% SubLarge
gammaClass(:,csubctr) =  ~ gammaClass(:,csub);  %% CTR
if ~isempty(evSup)

    %% SUPER thres
    temp = findSimEvent(cinj_ind,evSup,round(wind/1000/output.dt))'
    gammaClass(temp1(:,2),csup)= 1;  %% Sub
    temp = findSimEvent(cinj_indS,evSup,round(wind/1000/output.dt))'
    gammaClass(temp1(:,2),csupS)= 1;  %% SubSmall
    temp = findSimEvent(cinj_indL,evSup,round(wind/1000/output.dt))'
    gammaClass(temp1(:,2),csupL) = 1; %% SubLarge
    % Ctrl
    gammaClass(:,csupctr) = gammaClass(:,csup) & (~(gammaClass(:,csupS)&gammaClass(:,csupL))); %% CTR
end
%%% check it works till here.
sprintf('REPORT// GC: %d\tGCsub:%d\tGCctr:%d\tGCLarge: %d\tGCSmall: %d', length(gammaClass),sum(gammaClass(:,csub)),sum(gammaClass(:,csubctr)),sum(gammaClass(:,csubL)),sum(gammaClass(:,csubS) ))





%% find ind_x that occur in phase X, Y, Z to creat average EPSP during
%% those phases. (can do this by adding time to the output of find Sim
%% Event)

%%%%

% plot histogram of times of conductance relative to gamma
%** this is a control of to make sure data is uniform.

% plot cell voltage during large/small gamma
% plot cell voltage during conductance injection and large/small gamma
%% plot avg epsp during different phases.
sel_ind = find(gammaClass(:,csubctr));
emSUB_ctr = mean(events(sel_ind,:),1) ;
esSUB_ctr = std(events(sel_ind,:),1);
xtime = [1:length(emSUB_ctr)].*output.dt*1000;
figure(10)
clf
plot(xtime,emSUB_ctr)
% INTRACELLULAR CURRENT triggered on osc
tmSUB_ctr = mean(eventtrigger(sel_ind,:),1);
tsSUB_ctr = std(eventtrigger(sel_ind,:),1);
figure(11)
clf
% plot(xtime,tmSUB_ctr)
% plot(mean(eventtrigger(sel_ind,:)'))

plot(tmSUB_ctr,'-k','linewidth',2)
hold on


sel_ind = find(gammaClass(:,csubL));
tmSUB_L = mean(eventtrigger(sel_ind,:));
plot(tmSUB_L,'-r','linewidth',2)
sel_ind = find(gammaClass(:,csubS));
plot(mean(eventtrigger(sel_ind,:)),'-g','linewidth',2)


tmSUB_ctr = tmSUB_ctr - mean(tmSUB_ctr);
tmSUB_L = tmSUB_L - mean(tmSUB_L);
figure(12)
plot(xtime,tmSUB_ctr)
hold on
plot(xtime,tmSUB_L,'-k')

for i = 1:size(events,1)
  Vmax(i) = max(eventtrigger(i,340:540))    ;
end
s_Plot_by_Quad


% 
% emeanSUB_L = [mean(events(:,:),1) std(events(:,:),1)];
% emeanSUB_S = [mean(events(:,:),1) std(events(:,:),1)];
% 
% evcinjL
% evcinjS
% evctr
% 
% %% MAKE DIR
% fileindex =  output.sfilename(max(strfind(output.sfilename,'_'))+1:max(strfind(output.sfilename,'_'))+4);
% writedirheader = [writedirheader output.sfilename(1:max(strfind(output.sfilename,'_'))-1) '_' exptnum '\'];
% if ~isdir(writedirheader)
%     mkdir(writedirheader);
% end;
% 
% %% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figString = 'rawdata';
% figure(2);
% clf
% %% RAWDATA LFP
% subplot(2,1,1)
% plot([1:1/output.dt].*output.dt*1000,prodata(1:1/output.dt))
% %     plot(prodata(1:9999,1))
% hold on;
% plot (Strough(1:max(find(Strough < 1/output.dt ))).*output.dt*1000,thres*ones(1,[1:max(find(Strough < 1/output.dt))]),'.r')
% %     plot (timeinsw(1:max(find(trough < 9999 ))),thres*ones(1,[1:max(find(trough < 9999))]),'.r')
% axis tight
% ylim([min(prodata(.1*1/output.dt:1/output.dt)) max(prodata(.1*1/output.dt:1/output.dt))])
% xlim([100 1000]);
% %% RAWDATA INTRA
% subplot(2,1,2)
% %     plot(intradata(1:9999,1)*output.gain(1));
% plot([1:1/output.dt].*output.dt*1000,intradata(1:1/output.dt)+output.offset(1));
% xlabel('time (ms)')
% if bVC
%     ylabel('current (pA)')
% else
%     ylabel('voltage (mV)')
% end
% axis tight
% xlim([100 1000]);
% 
% if breport
%     dirtemp = 'REPORT';
%     figdesc = [figString fileindex];
%     savefigure(writedirheader,dirtemp,figdesc,savetype)
% end
% 
% 
% %% SAVE DATA
% odata = struct();
% odata.expttype = expttype;
% odata.exptnum = [exptnum '_' output.sfilename(1:strfind(output.sfilename,'.')-1)]; %% UNIQUE KEY
% odata.sfilename = output.sfilename;
% odata.dt = output.dt;
% odata.holdVoltage = v;
% odata.VC = bVC;
% odata.thres = thres;
% odata.LPFilter = LPFilter;
% odata.ISIrange = ISIrange;
% odata.analysisDate = datestr(now);
% 
% odata.prodata1sec = prodata(1:1/output.dt);
% odata.intradata1sec = intradata(1:1/output.dt)+output.offset(1);
% 
% 
