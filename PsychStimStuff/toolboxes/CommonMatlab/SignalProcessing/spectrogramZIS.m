function [P F T] = spectrogramZIS(S,dt,fmin,tres,bplot)
% function [P F T] = spectrogramZIS(S,dt,fmin,tres,bplot)
% Testingn ZIS ('zero interval subraction) time-frequency algorithm from
% Am J Physiol Regulatory Integrative Comp Physiol 291:1414-1429, 2006. First published Jul 6, 2006;
% Vitaliy Marchenko and Robert F. Rogers
% http://ajpregu.physiology.org/cgi/reprint/291/5/R1414?ijkey=


bdebug = 0; % in debug mode S is overwritten
if bdebug
    % zwin = 1/dt*.04
    % fmin = 5;%(Hz) minimum frequency resolution
    dt = 1e-4; % intersample interval
    % clear tres;

    t = [0:dt:.6];
    lent=length(t);
    S = [sin(2*pi*10*t(1:lent/3)) sin(2*pi*10*t(lent/3+1:lent*2/3))+sin(2*pi*80*t(lent/3+1:lent*2/3))  sin(2*pi*10*t(lent*2/3+1:end))+sin(2*pi*80*t(lent*2/3+1:end))+sin(2*pi*160*t(lent*2/3+1:end))];
    figure(10);
    plot(S)

end

if nargin<5
    bplot =0;
end

lenS = length(S) ;%
% determine the window length based on the above parameters
% NB the shorter the window relative tothe zeroed window the bettwer small changes in freq will be detected on the other hand the lower the resolution.
if isempty(fmin)
    fmin = 1/(lenS*dt); % default resolution is max resolution allowed by length ot data
end

if tres>=1/fmin
    error('Zeroed window is as large time window')
end
win = (fmin)/dt; %minimum length for resolution
if isempty(tres)
    zwin= win/5;
else
    zwin = tres/dt;
end
% windowing functions are not used because it is not clear how to implement
% them since they must be implemented the same way on both the zeroed and
% non zeroed time series.

% win = length(S)/5;
noverlap = []; nfft = [];

n = ceil(length(S)/zwin);
if ~isrowvector(S);    S = S';end
S = [S zeros(1,n*zwin-length(S)+win)]; % make S into an integer number of zwindows

mask = ones(1,win); mask(1:zwin)=0; clear P;clear T;
tic
for i=1:2*n-1
    seg = S(int32(1+round(zwin/2*(i-1))):int32(win +round(zwin/2*(i-1))));
    temp = seg.*mask;    
 [ttemp2 F] =  pwelch(seg,dpss(round(win),0.5),noverlap,nfft,1/dt);
 P(i,:)= ttemp2 -pwelch(temp,dpss(round(win),0.5),noverlap,nfft,1/dt);
     T(i)=(i-1)*zwin/2*dt; % time of the beginning of the spectrogram bin
     if bdebug
         subplot(2,1,1)
         plot(F,ttemp2); xlim([0 80])
         subplot(2,1,2)
         plot(F,P(i,:));         xlim([0 80])
         pause;
     end
        
end
%%cut it off after exceed win or move zwin in win
toc
[pf ff] = pwelch(seg,dpss(round(win),0.5),noverlap,nfft,1/dt)
%% plot
if bplot
figure(12);clf
ImageMatrix(T,F,P');
axis xy;  colorbar('off'); 
ylim([0 100])
end
%