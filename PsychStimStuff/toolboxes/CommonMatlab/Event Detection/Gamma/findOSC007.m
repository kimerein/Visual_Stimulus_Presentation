

% clear all
clear all
defineDir
warning off
expttype = 'KinateOsc';
writedirpath = [DATAANAL_DIR expttype '\CA3DY\'];
writedirheader = [DATAANAL_DIR expttype '\CA3DY\'];
breport = 0;
bsave = 0;
savetype = 'emf';
exptnum = 's2c2';
if ~exist('lastexptnum')
    lastexptnum = '';
end
bplot =1;
% exptnum = 's2c1';
% exptnum = '2s6c1';
% % xcolor = '-g'; v = -3; % IPSCs
xcolor = '-r'; v = -87; % EPSCs
% xcolor = '-c'; v = -62;% mixed
bVC = 1;
output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_10_22_0041ss.abf',-1,1);
% output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_09_21_0041_IN1LFP.abf',-1,1);
% output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_07_20_009s.abf',-1,1);

if ~strcmp(lastexptnum,exptnum)
    close all
end



thres = -0.01; %% MANUALLY DEFINE threshold for oscilation detection
% thres = -2; %% MANUALLY DEFINE threshold for oscilation detection
LPFilter = 60 ;%% low pass filter of LFP
ISIrange = [1/20 1/50]*1000; % ms
odata = struct();
odata.expttype = expttype;
odata.exptnum = [exptnum '_' output.sfilename(1:strfind(output.sfilename,'.')-1)]; %% UNIQUE KEY
odata.sfilename = output.sfilename;
odata.dt = output.dt;
odata.holdVoltage = v;
odata.VC = bVC;
odata.thres = thres;
odata.LPFilter = LPFilter;
odata.ISIrange = ISIrange;
odata.analysisDate = datestr(now);

fileindex =  output.sfilename(max(strfind(output.sfilename,'_'))+1:max(strfind(output.sfilename,'_'))+4);
writedirheader = [writedirheader output.sfilename(1:max(strfind(output.sfilename,'_'))-1) '_' exptnum '\'];
if ~isdir(writedirheader)
    mkdir(writedirheader);
end;
indexoffset = 1; b_intra= 0;


%%LFP
lfpdata = double(output.data(:,indexoffset+3:output.Nchan:end))*output.gain(2);
% lfpdata = lfpdata/mean(std(lfpdata));
lfpdata = (lfpdata-mean(mean(lfpdata)));
intradata = double(output.data(:,indexoffset+1:output.Nchan:end))*output.gain(1);
% intradata = intradata - mean(mean(intradata));
% intradata = intradata/min(std(intradata));

%% low pass filter LFP data
[B,A] = butter(2,2*LPFilter*output.dt,'low');

prodata =   filtfilt(B,A,lfpdata);
% prodata =   lfpdata;
bdebug = 0;
if bdebug
    plot(prodata(1:5000,2))
    hold on
    plot(t(1:5000,2),'r')
end
thres = -10;
%% find trough of each oscilation
%% reject events out side of ISIrange
%% extract WOI in Intradata after event
findOSC_detectOSC001
figString = 'rawdata';
figure(1);
clf
%% RAWDATA LFP
subplot(2,1,1)
plot([1:1/output.dt].*output.dt*1000,prodata(1:1/output.dt))
%     plot(prodata(1:9999,1))
hold on;
plot (Strough(1:max(find(Strough < 1/output.dt ))).*output.dt*1000,thres*ones(1,[1:max(find(Strough < 1/output.dt))]),'.r')
%     plot (timeinsw(1:max(find(trough < 9999 ))),thres*ones(1,[1:max(find(trough < 9999))]),'.r')
axis tight
ylim([min(prodata(.1*1/output.dt:1/output.dt)) max(prodata(.1*1/output.dt:1/output.dt))])
xlim([100 1000]);
%% RAWDATA INTRA
subplot(2,1,2)
%     plot(intradata(1:9999,1)*output.gain(1));
plot([1:1/output.dt].*output.dt*1000,intradata(1:1/output.dt)+output.offset(1));
xlabel('time (ms)')
if bVC
    ylabel('current (pA)')
else
    ylabel('voltage (mV)')
end
axis tight
xlim([100 1000]);

if breport
    dirtemp = 'REPORT';
    figdesc = [figString fileindex];
    savefigure(writedirheader,dirtemp,figdesc,savetype)
end

odata.prodata1sec = prodata(1:1/output.dt);
odata.intradata1sec = intradata(1:1/output.dt)+output.offset(1);


emean = mean(events(:,:),1);
estd = std(events(:,:),1);
xtime = [1:length(emean)].*output.dt*1000;
if bplot
    figure(10)
    clf
    orient tall
    subplot(2,1,1)
    figString = 'LFPdata';

    plot(ones(size(events,1),1)*xtime,events(:,:),'.k');
    hold on;
    plot(xtime,emean,'-b','LineWidth',2.5);
    plot(xtime,emean-estd,'--b','LineWidth',1.5);
    plot(xtime,emean+estd,'--b','LineWidth',1.5);
    axis tight
    set(gca,'Color','none','XGrid', 'on')
    %interosc interval
    subplot(2,1,2)
    [N,X]=hist(ISI,1000);
    hist(ISI,1000);
    hold on;
    line([ISIrange(1) ISIrange(2);ISIrange(1) ISIrange(2)],[max(N) max(N); 0 0])
    lim = axis;
    text(mean(lim(1:2)),lim(4)*.9,['Events:' num2str(Ncycles)]);
    % title('Histogram InterEvent Interval')
    xlim([0 200]);
    xlabel('time (ms)');
    set(gca,'Color','none','XGrid', 'on')
    if breport
        dirtemp = 'REPORT';
        figdesc = [figString fileindex];
        savefigure(writedirheader,dirtemp,figdesc,savetype)
    end

end
odata.emean = emean;

%  subplot(3,1,3)



% INTRACELLULAR CURRENT triggered on osc
tmean = mean(eventtrigger(:,:).*output.gain(1),1)+output.offset(1);
tstd = std(eventtrigger(:,:).*output.gain(1),1);
if bplot
    figure(11)
    figString = 'TrigIntra';
    orient tall
    clf
    plot(ones(size(eventtrigger,1),1)*xtime,eventtrigger(:,:).*output.gain(1)+output.offset(1),'.k');
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
odata.tmeanlagpk = (odata.tmeanlag - temp)*output.dt*1000; %% % PSP lag from peak of LFP (ms)
odata.tmeanlag = (odata.tmeanlag)*output.dt*1000-40; %% % PSP lag from alignment point (probably 20% of peak) (ms)

if bplot
    figure(12)
    figString = 'LFPIntra';
    orient tall
    clf
    plot(xtime,emean/norm(emean),'LineWidth',2.5);
    hold on
    plot(xtime,(temptmean-mean(temptmean))/max(temptmean),xcolor,'LineWidth',2.5)
    axis tight
    set(gca,'Color','none','XGrid', 'on')
    xlabel('time (ms)')
        title(['lag_{pk}:' num2str(odata.tmeanlagpk,'%1.1f')]);

    if breport
        dirtemp = 'REPORT';
        figdesc = [figString fileindex];
        savefigure(writedirheader,dirtemp,figdesc,savetype)
    end
end

figure(13)
figString = 'ExIn_LFPIntra';
plot(xtime,emean/max(emean),'LineWidth',2.5);
hold on
plot(xtime,(temptmean-mean(temptmean))/max(temptmean-mean(temptmean)),xcolor,'LineWidth',2.5)
plot(xtime,(temptmean-mean(temptmean))/max(temptmean-mean(temptmean)),xcolor,'LineWidth',2.5)

set(gca,'Color','none','XGrid', 'on');  orient tall
xlabel('time (ms)')
axis tight
title(['N: ' num2str(size(eventtrigger,1)) ' lag_{pk}:' num2str(odata.tmeanlagpk,'%1.1f')]);

if breport
    dirtemp = 'REPORT';
    figdesc = [figString];
    savefigure(writedirheader,dirtemp,figdesc,savetype)
end

% d = mIPSC/norm(mIPSC) + mEPSC/norm(mEPSC);
% plot(mIPSC/norm(mIPSC),'-c','LineWidth',2.5);
% plot(d,'k','LineWidth',2.5);

%
%
t =     xcov(prodata(:,2),30000,'coeff');
tt =     xcov(intradata(:,2),30000,'coeff');
if bplot
    figure(16)
    clf
    figString = 'AutoCor';

    plot(([1:size(t,1)]-((size(t,1)-1)/2 +1)).*output.dt.*1000,t)
    hold on
    plot(([1:size(tt,1)]-((size(tt,1)-1)/2 +1)).*output.dt.*1000,tt-.5,xcolor)
    axis tight
    set(gca,'Color','none','XGrid', 'on'); orient tall
    if breport
        dirtemp = 'REPORT';
        figdesc = [figString fileindex];
        savefigure(writedirheader,dirtemp,figdesc,savetype)
    end
end
odata.autocorrLFP = t;
odata.autocorrIntra = tt;

ttt =     xcov(intradata(:,2),prodata(:,2),10000,'coeff');
if bplot
    figure(17)
    figString = 'crosscorr';
    if ~strcmp(lastexptnum,exptnum)
        clf
    end
    set(gca,'Color','none','XGrid', 'on'); orient tall
    title1 = 'xcov LFP and Intracellular currents';
    hold on
    if strcmp(xcolor,'-g')
        plot(([1:size(ttt,1)]-((size(ttt,1)-1)/2 +1)).*output.dt.*1000,-ttt,xcolor)
    else
        plot(([1:size(ttt,1)]-((size(ttt,1)-1)/2 +1)).*output.dt.*1000,ttt,xcolor)
    end
    xlabel('ms');
    axis tight
%     if p < 0.05 then correlated
        if breport
            dirtemp = 'REPORT';
            figdesc = [figString];
            savefigure(writedirheader,dirtemp,figdesc,savetype)
        end
%     end

end
[r,p] = corrcoef(intradata(:,:),prodata(:,:));

odata.xcovLFPIntra = ttt;

%%% DATA SPECIFIC
% xn = reshape(intradata(:,1:end),1,[]);
% yn = reshape(prodata(:,1:end),1,[]); %% this is filtered
xn = reshape(intradata(:,2:end),1,[]);
yn = reshape(lfpdata(:,2:end),1,[]);

% %%% Power Spectrum
findOSC_s_powerspect
lowgamma= min(find(f>(1./ISIrange(1)*1000)));
highgamma = max(find(f<(1./ISIrange(2)*1000)));
odata.psdLFP = [f';Pyy]';
odata.psdIntra =[f';iPyy]';
odata.maxGPowLFP = max(Pyy(lowgamma:highgamma));
odata.maxGPowIntra = max(iPyy(lowgamma:highgamma));
odata.freqmaxGPowLFP = f(find(Pyy(lowgamma:highgamma)== odata.maxGPowLFP)+lowgamma);
odata.freqmaxGPowIntra = f(find(iPyy(lowgamma:highgamma)== odata.maxGPowIntra)+lowgamma);

%%% ADD PSD
odata


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
