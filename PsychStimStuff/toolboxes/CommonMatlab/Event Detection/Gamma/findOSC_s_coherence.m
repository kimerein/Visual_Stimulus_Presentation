nfft = 4*2048;
noverlap = 512;
nwin = hanning(1024);
fs = 1/output.dt;

[CxyCohere,F1] = mscohere(xn,yn,nwin,noverlap,nfft,fs);
% Compute the PSDs using Welch's method, which is similar to the
% MSCOHERE function
[Sxx,F] = pwelch(xn,nwin,noverlap,nfft,fs);
[Syy,F] = pwelch(yn,nwin,noverlap,nfft,fs);

% Convert one-sided PSDs to Power Spectrums
Sxx = [Sxx(1); (Sxx(2:end-1)/2); Sxx(end)].*fs;
Syy = [Syy(1); (Syy(2:end-1)/2); Syy(end)].*fs;

% Call CPSD function with similar arguments as the MSCOHERE
% function from above.
[Sxy,F] = cpsd(xn,yn,nwin,noverlap,nfft,fs);
Cxy = (abs(Sxy).^2)./(Sxx.*Syy);
phase = unwrap(angle(Sxy)*180/pi);
if bplot
    figure(18)
    figString = 'Coherence';
    if ~strcmp(lastexptnum,exptnum)
        clf
    end
    set(gca,'Position',[0 500 1000 1000],'Color','none','XGrid', 'on'); orient tall
    title('coherence LFP and Intracellular currents');
    % Parameters common to both CPSD and MSCOHERE:



    %%% PLOTTING COHERENCE
    subplot(1,2,1)
    set(gca,'Color','none','XGrid', 'on')
    % plot(f1,abs(a1),xcolor);
    % plot(1:length(Cxy), [CxyCohere Cxy]);% Compare the estimatessubplot(1,2,1)
    plot(F1, [CxyCohere],xcolor);% Compare the estimatessubplot(1,2,1)
    xlim([1 100]);
    xlabel('freq (Hz)');
    hold on;
    subplot(1,2,2)
    set(gca,'Color','none','XGrid', 'on')
    plot(F,phase,xcolor);
    xlim([1 100]) ;
    ylabel('phase (degrees)')
    xlabel('freq (Hz)')
    hold on;

    if breport
        dirtemp = 'REPORT';
        figdesc = [figString];
        savefigure(writedirheader,dirtemp,figdesc,savetype)
    end
end