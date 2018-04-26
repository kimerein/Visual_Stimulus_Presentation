figure(19)
set(gcf,'Position',[575 50 1000 800],'PaperUnits','normalized')
clf
subplot(2,1,1)
figString = 'PowerSpect';
% flim = 2048*4;
Hs=spectrum.periodogram;
msspectrum(Hs,xn,'Fs',1/output.dt)
hold on
Fs=spectrum.periodogram;
% set(Fs,'Color','r')
msspectrum(Fs,yn,'Fs',1/output.dt)
xlim([0 .1]);
set(gca,'Color','none');  
 subplot(2,1,2)
 %%Spectrogram
spectrogram(reshape(prodata,1,[]),1/output.dt/5,[],2^12,1/output.dt,'yaxis');
%  [S,F,T,P] =spectrogram(reshape(prodata,1,[]),1/output.dt/5,[],2^12,1/output.dt);
%  F = [0:2:1000];
%   [S,F,T,P] = spectrogram(reshape(prodata,1,[]),1/output.dt/5,[],F,1/output.dt)
%   title('Spectrogram')
   ylim([0 100])
caxis([min(caxis)/3.5 max(caxis)])
 colormap(jet)
     if breport
        dirtemp = 'REPORT';
        figdesc = [figString fileindex];
        savefigure(writedirheader,dirtemp,figdesc,savetype)
     end
    