savefilename = 'Z:\Data Analysis\mdata\O-2_121605_O01.mat';
%%SEEMS Correct .. but should be double checked with someone else
load(savefilename,...
    'Ts','processedFiles','Nspikes_in_File','extracted_spikes', 'EE_coord','-mat');
% l y is found by taking the 512-point fast Fourier transform (FFT): 
y = mean(extracted_spikes(:,:,3),2)';
Y = fft(y,1024);
%The power spectrum, a measurement of the power at various frequencies, is 
Pyy = Y.* conj(Y) / 1024;
% Graph the first 514 points (the other 512 points are redundant) on a meaningful frequency axis: 
f = Ts*(0:513)/1024;
figure;
plot(f,Pyy(1:514),'.r')
title('Frequency content of y')
xlabel('frequency (Hz)')
% 
% # N = 2^9;     % number of points to analyze
% # c = fft(y(start:start+N-1))/N     % compute fft of sound data, take a look at the output.
% # p = 2*abs(c(2:N/2));     % compute power at each frequency, see the note at the end.
% # f = (1:N/2-1)*Fs/N;     % frequency corresponding to p
% 
% # semilogy(f,p)     % what is semilogy? 

%%%%%%%%%%% OR

    Hs=spectrum.periodogram;
%     Hs = spectrum.welch;
% psd(Hs,prodata(:,2),'Fs',1/output.dt)
msspectrum(Hs,prodata(:,2),'Fs',1/output.dt)

%%%%%%%%%%%%%
    %% FFT over time
%     figure(25);
%      spectrogram(prodata(:,2),512,1,flim,1/output.dt)
%      xlim([0 500])