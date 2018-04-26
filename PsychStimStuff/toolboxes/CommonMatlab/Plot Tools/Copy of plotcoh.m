function [CxyCohere F1 phase F] = plotcoh(x,y,dt,fnum,prange,bplotphase)
% function [CxyCohere F1 phase F] = plotcoh(x,y,dt,fnum,prange,bplotphase)
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
% NOTE: window and overlap parameters have been tweek for looking at gamma
% osc.  But I don't really understand how to properly optimize them.  
% Different  nfft and window size can have large effects on the phase plot
% in particular
%% defaults
if nargin < 4 || isempty(fnum)
fnum = 3000;
end
if nargin <5 || isempty(prange)
prange = [0 100];
end
if nargin <6
bplotphase = 1;
end
%% parameters
nfft = 8*2048;
nfftp = 4*2048;
noverlap = 2^8;
% if 1/dt > length(x) |1/dt > length(y)
%     nwin = hamming(2^nextpow2(min(length(x),length(y))));
% else
%     nwin = hamming(2^nextpow2(1/dt)); % choose resolution of 1Hz (unless x or y are too short)
% end
nwin = hamming(2^12);
nwinp = hamming(2^10); % 
 %% NOTE: for reasons I don't understand phase gets very noisy when nwin is
 %% large. But Coh has better resolution.
 %% Need to undeterstand better how to chose the window esspecially with
 %% regard to extracting phase.
fs = 1/dt;
%% Get amp
d = find(size(x) ==min(size(x))); %get the right dimension
if d ~=2
    x = x';
    y = y';
end
for i = 1:size(x,2)
    [CxyCohere(:,i),F1] = mscohere(x(:,i),y(:,i),nwin,noverlap,nfft,fs);

    %% Get Phase
    % Compute the PSDs using Welch's method, which is similar to the
    % MSCOHERE function
    [Sxx(:,i) ,F] = pwelch(x(:,i),nwinp,noverlap,nfftp,fs);
    [Syy(:,i) ,F] = pwelch(y(:,i),nwinp,noverlap,nfftp,fs);
    % Call CPSD function with similar arguments as the MSCOHERE
    % function from above.
    [Sxy(:,i) ,F] = cpsd(x(:,i),y(:,i),nwinp,noverlap,nfftp,fs);
end
if size(x,2) >1
    CxyCohere = mean(CxyCohere,2);
    Sxx = mean(Sxx,2);
    Syy = mean(Syy,2);
    Sxy = mean(Sxy,2);
end
% Convert one-sided PSDs to Power Spectrums
Sxx = [Sxx(1); (Sxx(2:end-1)/2); Sxx(end)].*fs;
Syy = [Syy(1); (Syy(2:end-1)/2); Syy(end)].*fs;

% Cxy = (abs(Sxy).^2)./(Sxx.*Syy); %% to manually get amplitude
phase = unwrap(angle(Sxy)*180/pi);

%% Plot
figure(fnum);
% set(gca,'Position',[0 500 1000 1000],'Color','none','XGrid', 'on');
% subplot(1,2,1)
plot(F1, [CxyCohere]);
title('Coherence');
%     plot(F1, Cxy); %% can alternatively use manually computed value
%     (similar but not identical)
xlim(prange);
xlabel('freq (Hz)');
hold on;
plotset(1);
%
if bplotphase
    figure(fnum+1);
    % set(gca,'Position',[0 1500 1000 1000],'Color','none','XGrid', 'on');
    % subplot(1,2,2)
    plot(F,phase);
    xlim(prange) ;
    ylabel('phase (degrees)')
    xlabel('freq (Hz)')
    hold on;
    plotset(1);
end