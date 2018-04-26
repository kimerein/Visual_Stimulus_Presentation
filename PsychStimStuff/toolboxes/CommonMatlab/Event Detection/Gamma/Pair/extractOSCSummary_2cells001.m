
defineDir
expttype = 'KinateOsc';
writedirpath = [DATAANAL_DIR expttype '\'];
writedirpath = [DATAANAL_DIR expttype '\'];
writedirheader = [DATAANAL_DIR expttype '\'];
breport = 1;
bsave = 1;
savetype = 'emf';

%% USER DEFINE
exptnum = 's4';
indexoffset = 1;
xcolor = '-k'; 
bVC = 1;


% output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_07_20_0011.abf',-1,1);
output = importStruct_abf('C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2006_09_18_0032s.abf',-1,1);

%%% names do not reflect reality .. often channel 1 and 2 are both
%%% intracellular
lfpdata = double(output.data(:,indexoffset+2:output.Nchan:end))*output.gain(2);
intradata = double(output.data(:,indexoffset+1:output.Nchan:end))*output.gain(1);


figure(30)
clf
set(gcf,'Position',[580 50 1000 700],'PaperUnits','normalized');    hold on;    orient tall
subplot(2,3,[1:3])
plot([1:1/output.dt].*output.dt*1000,lfpdata(1:1/output.dt)-mean(lfpdata(.1*1/output.dt:1/output.dt)),'r','LineWidth',2)
hold on;
plot([1:1/output.dt].*output.dt*1000,intradata(1:1/output.dt)-mean(intradata(.1*1/output.dt:1/output.dt)),'b','LineWidth',2)
axis tight
% ylim([min(lfpdata(.1*1/output.dt:1/output.dt)) max(lfpdata(.1*1/output.dt:1/output.dt))])
xlim([100 1000]);
title(output.sfilename,'Interpreter','none')

clear plfp; clear pintra;
for i=1:(min(size(lfpdata,2),40))
    [plfp(i,:) f] = PWELCH(lfpdata(:,i)-mean(lfpdata(:,i)),[],[],[],1/output.dt);
    [pintra(i,:) f] = PWELCH((intradata(:,i)-mean(intradata(:,i))),[],[],[],1/output.dt);
%       [pintra(i,:) f] = PWELCH((intradata(:,i)-mean(intradata(:,i)))/std(intradata(:,1)),[],[],[],1/output.dt);
      
      %% don't really know why I have to do this... seems to explode if
      
      %% it is not scaled down
end
% Pyy = 10.^mean(plfp,1)-1;
% iPyy = 10.^mean(pintra,1)-1;
Pyy = mean(plfp,1);
iPyy = mean(pintra,1);
clear plfp; clear pintra;
subplot(2,3,4)
plot(f,Pyy,'r','LineWidth',2) %%LFP
hold on
plot(f,iPyy,'b','LineWidth',2)
axis tight
xlim([0 100])
% ylim([0 max(Pyy(1:100))]*1.05);
xlabel('frequency (Hz)')
ylabel('pA^2/Hz')
title('PSD')

% ttt =     xcorr(intradata(:,2),lfpdata(:,2),1000,'unbiased');
ttt =     xcov(intradata(:,2),lfpdata(:,2),10000,'coeff');
subplot(2,3,5)
x = ([1:size(ttt,1)]-((size(ttt,1)-1)/2 +1)).*output.dt.*1000;
plot(x,ttt,xcolor,'LineWidth',2)
hold on
axis tight
xlim([-100 100])
temp = ttt((size(ttt,1)-1))/abs(ttt((size(ttt,1)-1)));
title(['CrossCov' ' peak: ' num2str(max(abs(ttt))*temp,'%1.2f') ' lag: ' num2str(x(find(ttt==(max(abs(ttt))*temp))),'%1.2f') 'ms']);


subplot(2,3,6)
xn = reshape(intradata(:,2:end),1,[]);
yn = reshape(lfpdata(:,2:end),1,[]);
nfft = 4*2048;
noverlap = 512;
nwin = hanning(1024);
fs = 1/output.dt;
[CxyCohere,F1] = mscohere(xn,yn,nwin,noverlap,nfft,fs);
plot(F1, [CxyCohere],xcolor,'LineWidth',2);% Compare the estimatessubplot(1,2,1)
xlim([1 100]);
xlabel('freq (Hz)');
hold on
axis tight
xlim([0 100])
title('Coherence')

figString = ['Summary_' exptnum '_' output.sfilename(1:max(strfind(output.sfilename,'.'))-1)]
if breport
    dirtemp = 'Pairs';
    figdesc = [figString];
    savefigure(writedirheader,dirtemp,figdesc,savetype)
end
% subplot(2,2,4)
% % Compute the PSDs using Welch's method, which is similar to the
% % MSCOHERE function
% [Sxx,F] = pwelch(xn,nwin,noverlap,nfft,fs);
% [Syy,F] = pwelch(yn,nwin,noverlap,nfft,fs);
% 
% % Convert one-sided PSDs to Power Spectrums
% Sxx = [Sxx(1); (Sxx(2:end-1)/2); Sxx(end)].*fs;
% Syy = [Syy(1); (Syy(2:end-1)/2); Syy(end)].*fs;
% 
% % Call CPSD function with similar arguments as the MSCOHERE 
% % function from above.
% [Sxy,F] = cpsd(xn,yn,nwin,noverlap,nfft,fs);
% Cxy = (abs(Sxy).^2)./(Sxx.*Syy);
% phase = unwrap(angle(Sxy)*180/pi);
% set(gca,'Color','none','XGrid', 'on')
% plot(F,phase,xcolor);
% xlim([1 100]) ;
% ylabel('phase (degrees)')
% xlabel('freq (Hz)')
% hold on;

