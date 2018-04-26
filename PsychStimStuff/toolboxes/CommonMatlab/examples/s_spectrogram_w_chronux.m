% example spectrogram with chronux
movingwin = [0.2 0.1]; %window, windowstep
sparams.Fs = 1/dt;
[S,t,f]=mtspecgramc(df(:,:)',movingwin,sparams);

[temp t] = convert_multisweep_spectrum(S,t);

figure(1);clf;
% subplot(2,1,1)
plotmatrix(temp,t,f)
ylim([0 100])
