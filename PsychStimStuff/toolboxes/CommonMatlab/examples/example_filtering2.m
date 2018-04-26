% Example FFT  - > filter IFFT
% and band pass filtering

width = size(pd,1);
g=fft(pd,width)

mag_g=abs(g);   
phs_g=angle(g);
x = [1:width];
m = 31; v=6;
gaus = (1/(v*sqrt(2*3.142))).*exp(-(x-m).^2./(2*v^2));
gaus = 1 - gaus./max(gaus);

figure(2);
plot(mag_g,'b');
mag_g = gaus.*mag_g';
plot(mag_g,'r');
%gauss=p
%gauss=gauss/max(gauss)
%window=ones(size(fft))-gauss
%mag*window

modified=mag_g.*(cos(phs_g')+(i.*sin(phs_g')));  %these lines of code
%are the same
%mod=magim.*exp(i*newPhase);

filtdata = ifft(modified, size(pd,1),'symmetric');

% mag_g*exp(-pi*phase)
figure(1);
hold off
plot(xtime,pd);
hold on;
plot(xtime,filtdata,'r');

%%% try filtering with filter
% % 
% high = 2*50*samplerate;
% low = 2*20*samplerate;
% 
% % [B,A] = butter(3,high,'high');
% % [B,A] = butter(3,[low high],'stop');
% [B,A] = ellip(3,1,100,[low high],'stop');
% L1 = filtfilt(B,A,double(pd));
% % 
% % low = 2*50*samplerate;
% % [B,A] = butter(2,low,'low');
% % L2 = filtfilt(B,A,double(L1));
% 

