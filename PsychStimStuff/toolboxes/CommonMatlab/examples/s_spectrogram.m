% Fdata.exptname = ['V100107c1'];
Fdata.sName = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\2007_09_30_0009.abf';
Fdata.sChan = {'IN 1      ','IN 15     '};
Fdata.sweeps = [];
Fdata.start = 7*60;
Fdata.stop = 10*60;

%%
[d dt analyzedTime] =loadEPfiles(Fdata);
% do subsampling here.
d = d(1:10:end,:); dt = dt*10;
%%
nwin = round(400e-3/dt); % time window to analyze
[win,v] = dpss(nwin,.5); % create slepian window
Res = 2; nfft = 2^nextpow2((1/dt/2)/Res);   % define resolution of fft 
% (note can't be higher than allowed by the length of the win)
[ S F T P] = spectrogram(d,win,[],nfft,1/dt,'yaxis');
figure;
a = max(find(F<= 5)); b = max(find(F<= 120)); % plot only subset of frequencies
surf(T,F(a:b),10*log10(abs(P(a:b,:))),'EdgeColor','none');
axis xy; axis tight; colormap(jet); view(0,90);
xlabel('Time'); ylabel('Frequency (Hz)');
title('Spectrogram')

%% coherence gram
figure;
[P F T] = mtchg(d,nfft,1/dt,nwin,[],1,1);% coherence gram
a = max(find(F<= 5)); b = max(find(F<= 120)); % plot only subset of frequencies
[H X Y ] = ImageMatrix(T,F(a:b),10*log10(P(a:b,:,:,:)));
for i = 1: length(H); subplot(H(end-i+1)); axis xy; view(0,90); colorbar('off'); end
xlabel('Time'); ylabel('Frequency (Hz)');
title('Spectrogram/Coherogram')

