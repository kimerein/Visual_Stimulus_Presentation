% clear all
defineDir
warning off
expttype = 'KinateOsc';
writedirpath = [DATAANAL_DIR expttype '\'];
writedirheader = [DATAANAL_DIR expttype '\'];
breport = 1;
bsave = 0;
savetype = 'emf'
exptnum = 's2c1';
% output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_07_20_0011.abf',-1,1);
output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_07_20_009s.abf',-1,1);


% xcolor = '-g'; v = 7.5; % IPSCs
% xcolor = '-r'; v = -87; % EPSCs
xcolor = '-c'; v = -48;% mixed
bVC = 1;

thres = -.3; %% MANUALLY DEFINE threshold for oscilation detection
 LPFilter = 60 %% low pass filter of LFP
odata = struct();
odata.expttype = expttype;
odata.exptnum = [exptnum '_' output.sfilename(1:strfind(output.sfilename,'.')-1)]; %% UNIQUE KEY
odata.sfilename = output.sfilename;
odata.dt = output.dt;
odata.holdVoltage = v;
odata.VC = bVC;
odata.thres = thres;
odata.LPFilter = LPFilter;
odata.analysisDate = datestr(now);

fileindex =  output.sfilename(max(strfind(output.sfilename,'_'))+1:max(strfind(output.sfilename,'_'))+4);
writedirheader = [writedirheader output.sfilename(1:max(strfind(output.sfilename,'_'))-1) '_' exptnum '\'];
if ~isdir(writedirheader)
    mkdir(writedirheader);
end;
indexoffset = 1; b_intra= 0;

  
%%LFP
lfpdata = double(output.data(:,indexoffset+2:output.Nchan:end));
lfpdata = lfpdata/min(std(lfpdata));
lfpdata = lfpdata-mean(mean(lfpdata));
intradata = double(output.data(:,indexoffset+1:output.Nchan:end));
intradata = intradata - mean(mean(intradata));
% intradata = intradata/min(std(intradata));

%% low pass filter LFP data
      [B,A] = butter(2,2*LPFilter*output.dt,'low');
prodata =   filtfilt(B,A,lfpdata);
bdebug = 0;
if bdebug
 plot(lfpdata(1:5000,2))
 hold on
 plot(t(1:5000,2),'r')
end
%% find trough of each oscilation
  %% threshold
a= prodata < thres;
a = [diff(prodata); zeros(1,size(prodata,2))].*a;
j = 1;
    for (i = 2: size(a,1)*size(a,2))
        if a(i) > 0 & a(i-1) <= 0 % if trough
            if min(abs(a(i)),abs(a(i-1)))==abs(a(i))
                trough(j) = i;
                baug = 1;
            else
                trough(j) = i-1;
                baug = 1;
            end
            if baug
                j = j+1;
                %% % plot to debug detection of troughs
                debug = 0;
                if debug
                    figure(100);
                    plot(a(1:19999),'k') %%% DATA SPECIFIC
                    hold on
                    plot(prodata(1:19999),'b')
                    plot(trough(j),.5,'.r')
                    pause;
                end
            end
        end
    end


%% Extract event
WOI = (1/output.dt/1000).*[40 10];
events = zeros(size(trough,1),sum(WOI)+1); eventtrigger = zeros(size(trough,1),sum(WOI)+1); 
for i = 1: size(trough,2)
    try
        %         check same sweep
        if int32((trough(i) - WOI(1))/size(prodata,1)) == int32((trough(i) + WOI(2))/size(prodata,1))
            eventsw = floor((trough(i) - WOI(1))/size(prodata,1));
            timeinsw(i) = trough(i) - (eventsw)*size(prodata,1);
            events(i,:) = prodata(trough(i)-WOI(1):trough(i)+WOI(2));
            eventtrigger(i,:) = intradata(trough(i)-WOI(1):trough(i)+WOI(2));
        end
        %         plot(events(i,:))
        %         pause;

    catch
        i
    end
end

figString = 'rawdata';
figure(1);
clf

%% RAWDATA LFP
subplot(2,1,1)
plot(prodata(1:1/output.dt))
%     plot(prodata(1:9999,1))
hold on;
plot (trough(1:max(find(trough < 1/output.dt ))),thres*ones(1,[1:max(find(trough < 1/output.dt))]),'.r')
%     plot (timeinsw(1:max(find(trough < 9999 ))),thres*ones(1,[1:max(find(trough < 9999))]),'.r')
axis tight
%% RAWDATA INTRA
subplot(2,1,2)
%     plot(intradata(1:9999,1)*output.gain(1));
plot([1:1/output.dt].*output.dt*1000,intradata(1:1/output.dt)*output.gain(1)+output.offset(1));
xlabel('time (ms)')
if bVC
    ylabel('current (pA)')
else
    ylabel('voltage (mV)')
end
axis tight
if breport
    dirtemp = 'REPORT';
    figdesc = [figString fileindex];
    savefigure(writedirheader,dirtemp,figdesc,savetype)
end

odata.prodata1sec = prodata(1:1/output.dt);
odata.intradata1sec = intradata(1:1/output.dt)*output.gain(1)+output.offset(1);


figure(10)
clf
orient tall
subplot(2,1,1)
figString = 'LFPdata';
emean = mean(events(:,:),1);
estd = std(events(:,:),1);
xtime = [1:length(emean)].*output.dt*1000;
plot(ones(size(events,1),1)*xtime,events(:,:),'.k');
hold on;
plot(xtime,emean,'-b','LineWidth',2.5);
plot(xtime,emean-estd,'--b','LineWidth',1.5);
plot(xtime,emean+estd,'--b','LineWidth',1.5);
axis tight
set(gca,'Color','none','XGrid', 'on')
%interosc interval
subplot(2,1,2)
aa = diff(trough)*output.dt*1000;
hist(aa,1000);
% title('Histogram InterEvent Interval')
xlim([0 200]);
xlabel('time (ms)');
set(gca,'Color','none','XGrid', 'on')
    if breport
        dirtemp = 'REPORT';
        figdesc = [figString fileindex];
        savefigure(writedirheader,dirtemp,figdesc,savetype)
    end

    odata.emean = emean;
    
%  subplot(3,1,3)



% INTRACELLULAR CURRENT triggered on osc
figure(11)
figString = 'TrigIntra';
orient tall
tmean = mean(eventtrigger(:,:).*output.gain(1),1);
tstd = std(eventtrigger(:,:).*output.gain(1),1);
clf
plot(ones(size(eventtrigger,1),1)*xtime,eventtrigger(:,:).*output.gain(1),'.k');
hold on;
plot(xtime,tmean,xcolor,'LineWidth',2.5);
plot(xtime,tmean-tstd,'--r','LineWidth',1.5);
plot(xtime,tmean+tstd,'--r','LineWidth',1.5);
axis tight
set(gca,'Color','none','XGrid', 'on')
xlabel('time (ms)')
if bVC
    ylabel('current (pA)')
else
    ylabel('voltage (mV)')
end
if breport
    dirtemp = 'REPORT';
    figdesc = [figString fileindex];
    savefigure(writedirheader,dirtemp,figdesc,savetype)
end

odata.tmean = tmean;
    
    %% flip in IPSC
    temptmean = tmean;
if strcmp(xcolor,'-g')
    temptmean = -1*temptmean;
end
%%% COLLECT SUMMARY DATA
   odata.tmeanmin = min(temptmean ); %% amplitude of PSP minimum
   odata.tmeanstd = tstd(find(temptmean==odata.tmeanmin)); %%std at PSP minimum
   odata.tmeanlag = find(temptmean==odata.tmeanmin);
    temp = find(emean == min(emean))
    odata.tmeanlag = (odata.tmeanlag - temp)*output.dt*1000; %% PSP Lag (ms)
    
figure(12)
figString = 'LFPIntra';
orient tall
clf
plot(xtime,emean/norm(emean),'LineWidth',2.5);
hold on
plot(xtime,tmean/norm(tmean),xcolor,'LineWidth',2.5)
axis tight
set(gca,'Color','none','XGrid', 'on')
xlabel('time (ms)')
    if breport
        dirtemp = 'REPORT';
        figdesc = [figString fileindex];
        savefigure(writedirheader,dirtemp,figdesc,savetype)
    end
% title(['lag:' num2str(tmeanlag,'%f1.1')]);


figure(13)
if ~strcmp(lastexptnum,exptnum)
clf
end
figString = 'ExIn_LFPIntra';
plot(xtime,emean/norm(emean),'LineWidth',2.5);
hold on
plot(xtime,temptmean/norm(temptmean),xcolor,'LineWidth',2.5)
set(gca,'Color','none','XGrid', 'on');  orient tall
xlabel('time (ms)')
    if breport
        dirtemp = 'REPORT';
        figdesc = [figString];
        savefigure(writedirheader,dirtemp,figdesc,savetype)
    end

% d = mIPSC/norm(mIPSC) + mEPSC/norm(mEPSC);
% plot(mIPSC/norm(mIPSC),'-c','LineWidth',2.5);
% plot(d,'k','LineWidth',2.5);
axis tight


     
figure(16)
clf
figString = 'AutoCor';
t =     xcorr(prodata(:,2),10000,'coeff');
plot(t)
hold on
tt =     xcorr(intradata(:,2),10000,'coeff');
plot(tt-.5,xcolor)
axis tight
set(gca,'Color','none','XGrid', 'on'); orient tall
if breport
    dirtemp = 'REPORT';
    figdesc = [figString fileindex];
    savefigure(writedirheader,dirtemp,figdesc,savetype)
end
odata.autocorrLFP = t;
odata.autocorrIntra = tt;

    
 figure(17)
figString = 'crosscorr';
if ~strcmp(lastexptnum,exptnum)
clf
end
set(gca,'Color','none','XGrid', 'on'); orient tall
title = 'xcorr LFP and Intracellular currents';
hold on
ttt =     xcorr(intradata(:,2),prodata(:,2),10000,'coeff');
if strcmp(xcolor,'-g')
plot(-ttt,xcolor)
else
    plot(ttt,xcolor)
end
axis tight
%  [r,p] = corrcoef(intradata(:,:),prodata(:,:));
% if p < 0.05 then correlated
if breport
    dirtemp = 'REPORT';
    figdesc = [figString];
    savefigure(writedirheader,dirtemp,figdesc,savetype)
end
odata.xcorrLFPIntra = ttt;

%%% DATA SPECIFIC
% xn = reshape(intradata(:,1:end),1,[]);
% yn = reshape(prodata(:,1:end),1,[]);
xn = reshape(intradata(:,2:40),1,[]);
yn = reshape(prodata(:,2:40),1,[]);

% %%% Power Spectrum
findOSC_s_powerspect
lowgamma= min(find(f>20));
highgamma = max(find(f<50));
odata.psdLFP = Pyy;
odata.psdIntra = iPyy;
odata.maxGPowLFP = max(Pyy(lowgamma:highgamma));
odata.maxGPowIntra = max(iPyy(lowgamma:highgamma));
odata.freqmaxGPowLFP = f(find(Pyy(lowgamma:highgamma)== odata.maxGPowLFP)+lowgamma);
odata.freqmaxGPowIntra = f(find(iPyy(lowgamma:highgamma)== odata.maxGPowIntra)+lowgamma);

     %%% ADD PSD 
%      odata


% %%Coherence
%% Manual Coherence copied from "Spectral Estimation Method" in Matlab Doc
%% cpsd used to get phase
findOSC_s_coherence
%% COLLECT DATA for summary
odata.cohpeak = max(CxyCohere); %% value of max coherence
odata.freqcohpeak = F1(find(CxyCohere==odata.cohpeak));  %% Peak coherence frequency
odata.freqfwhm = fwhm2(F1,CxyCohere); %% FWHM of peak
odata.phasecohpeak = phase(find(F==odata.freqcohpeak)); %% phaselag at peak
odata.cohlag =  odata.phasecohpeak/180* 1/odata.freqcohpeak*1000;%% lag in ms
odata.cohLFPIntra = [F1 CxyCohere];
odata.phaseLFPIntra = [F phase];
%% AUTO coherence .. not sure this makes sense
% figure(19)
% [a f]= mscohere(intradata(:,2),intradata(:,3),hanning(1024),512,2048*4,1/output.dt);
% plot(f,a,'-r');
% hold on
% xlim([1 100])
if bsave
    sumdata = 'SummaryData.mat';
    if ~isempty(dir([writedirpath sumdata]))
        load([writedirpath sumdata]);
        dupfound = 0;
        for i=1:size(summarydata,2)  %% search for duplicates of the same analysis
            if strcmp(summarydata(i).exptnum,odata.exptnum)
                dupfound = 1;
                break;
            end
        end
        if dupfound
            summarydata(i) = odata;
        else
            ii= 1+ size(summarydata,2)
            summarydata(ii) =  odata;
        end
    else
        summarydata = odata;
    end


    save([writedirpath sumdata],'summarydata');
end
lastexptnum = exptnum;