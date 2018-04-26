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
%% Problem here is that I am not sure what the units are.
%%Pwelch outputs dB/Hz. prodata should be in mV but is in arbitary units  because importStruct_abf somehow screws up the scaling.
%% So prodata is in units of std
if bplot
figure(15)
clf
figString = 'PowerSpect';
set(gcf,'Position',[575 50 1000 800],'PaperUnits','normalized')
subplot(2,2,1)
plot(f,Pyy) %%LFP
xlim([0 100])
ylim([0 max(Pyy(1:100))]*1.05);
xlabel('frequency (Hz)')
ylabel('mV^2/Hz')
subplot(2,2,3)

plot(f,iPyy,xcolor)
xlim([0 100])
% ylim([0 max(iPyy(1:100))]*1.05/norm(iPyy));
ylim([0 max(iPyy(1:100))]*1.05);
xlabel('frequency (Hz)')
ylabel('pA^2/Hz')

subplot(2,2,2)
set(gca,'Color','none','XGrid', 'on')
 %%Spectrogram
spectrogram(reshape(lfpdata,1,[]),1/output.dt/5,[],2^12,1/output.dt,'yaxis');
%  [S,F,T,P] =spectrogram(reshape(prodata,1,[]),1/output.dt/5,[],2^12,1/output.dt);
%  F = [0:2:1000];
%   [S,F,T,P] = spectrogram(reshape(prodata,1,[]),1/output.dt/5,[],F,1/output.dt)
%   title('Spectrogram')
   ylim([0 100])
caxis([(max(caxis)-(max(caxis)-min(caxis))/5) max(caxis)])
subplot(2,2,4)
set(gca,'Color','none','XGrid', 'on')
 %%Spectrogram INTRA
spectrogram(reshape(intradata,1,[]),1/output.dt/5,[],2^12,1/output.dt,'yaxis');
%  [S,F,T,P] =spectrogram(reshape(prodata,1,[]),1/output.dt/5,[],2^12,1/output.dt);
%  F = [0:2:1000];
%   [S,F,T,P] = spectrogram(reshape(prodata,1,[]),1/output.dt/5,[],F,1/output.dt)
%   title('Spectrogram')
   ylim([0 100])
caxis([(max(caxis)-(max(caxis)-min(caxis))/5) max(caxis)])
    if breport
        dirtemp = 'REPORT';
        figdesc = [figString fileindex];
        savefigure(writedirheader,dirtemp,figdesc,savetype);
    end
end