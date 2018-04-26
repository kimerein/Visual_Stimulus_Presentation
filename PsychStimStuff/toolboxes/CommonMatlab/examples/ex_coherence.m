function [CxyCohere phase F] = plotcoh(x,y,dt,fnum,prange)
% function [CxyCohere phase F] = plotcoh(x,y,dt,fnum,prange)
% DESC:
% runs mscohere and manually calculates phase: 
%   Cxy = (abs(Sxy).^2)./(Sxx.*Syy);
% phase = unwrap(angle(Sxy)*180/pi);
%
% INPUT
%   x, y:    vectors
%   dt:      sample time (s)
%   fnum:    number of figure to create
%   prange:  range to plot (ms) default [-100 100]
% OUTPUT
%   CxyCohere:    real value
%   phase:   lagtimes
%   F:     frequency
%
% BA012707
%% defaults
fnum = 100;
prange = [0 100]
%% parameters
nfft = [];%4*2048;
noverlap = [];%512;
if 1/dt > length(x) |1/dt > length(y)
    nwin = min(length(x),length(y));
else
    nwin = min(1/dt); % choose resolution of 1Hz (unless x or y are too short)
end
fs = 1/dt;
%% Get amp
[CxyCohere,F1] = mscohere(xn,yn,nwin,noverlap,nfft,fs);

%% Get Phase
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
% Cxy = (abs(Sxy).^2)./(Sxx.*Syy); %% to manually get amplitude
phase = unwrap(angle(Sxy)*180/pi);

%% Plot
    figure(fnum)
%     figString = 'Coherence';
    set(gca,'Position',[0 500 1000 1000],'Color','none','XGrid', 'on'); 
    title('Coherence');
    % Parameters common to both CPSD and MSCOHERE:
    %%% PLOTTING COHERENCE
    subplot(1,2,1)
    set(gca,'Color','none','XGrid', 'on')
    plot(F1, [CxyCohere]);% Compare the estimatessubplot(1,2,1)
%     plot(F1, Cxy);
    xlim(prange);
    xlabel('freq (Hz)');
    hold on;
    subplot(1,2,2)
    set(gca,'Color','none','XGrid', 'on')
    plot(F,phase);
    xlim(prange) ;
    ylabel('phase (degrees)')
    xlabel('freq (Hz)')
    hold on;
