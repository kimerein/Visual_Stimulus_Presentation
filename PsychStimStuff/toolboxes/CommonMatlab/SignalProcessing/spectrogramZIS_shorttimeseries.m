function [P F T] = spectrogramZIS(S,dt,fmin,tres,bplot)
% function [P F T] = spectrogramZIS(S,dt,fmin,tres,bplot)
% Testingn ZIS ('zero interval subraction) time-frequency algorithm from
% Am J Physiol Regulatory Integrative Comp Physiol 291:1414-1429, 2006. First published Jul 6, 2006;
% Vitaliy Marchenko and Robert F. Rogers
% http://ajpregu.physiology.org/cgi/reprint/291/5/R1414?ijkey=

% fmin = 5;%(Hz) minimum frequency resolution
% dt = 1e-4; % intersample interval
% clear tres;

bdebug = 0; % in debug mode S is overwritten
if nargin<5
    bplot =0;
end

lenS = length(S) ;%
% determine the window length based on the above parameters
% NB the shorter the window relative tothe zeroed window the bettwer small changes in freq will be detected on the other hand the lower the resolution.
if isempty(fmin)
    fmin = 1/(lenS*dt); % default resolution is max resolution allowed by length ot data
end
if isempty(tres)
    tres= 1/fmin/5;
end

if tres>=1/fmin
    error('Zeroed window is as large time window')
end
win = 1/(fmin)/dt; %minimum length for resolution
zwin = tres/dt;
if bdebug
    % zwin = 1/dt*.04
    t = [0:dt:.6];
    lent=length(t);
    S = [sin(2*pi*10*t(1:lent/3)) sin(2*pi*10*t(lent/3+1:lent*2/3))+sin(2*pi*80*t(lent/3+1:lent*2/3))  sin(2*pi*10*t(lent*2/3+1:end))+sin(2*pi*80*t(lent*2/3+1:end))+sin(2*pi*160*t(lent*2/3+1:end))];
    figure(10);
    plot(S)

end
% windowing functions are not used because it is not clear how to implement
% them since they must be implemented the same way on both the zeroed and
% non zeroed time series.

% win = length(S)/5;
noverlap = []; nfft = [];
[p0 F] = pwelch(S,round(win),noverlap,nfft,1/dt);
n = ceil(length(S)/zwin);
if ~isrowvector(S);    S = S';end
S = [S zeros(1,n*zwin-length(S))]; % make S into an integer number of zwindows

tic
for i=1:2*n-1
    temp = S;
    temp(1+round(zwin/2)*(i-1):zwin+round(zwin/2)*(i-1))=0;
    p(i,:) = pwelch(temp,round(win),noverlap,nfft,1/dt);
%     P(i,:) = 
     T(i)=(i-1)*zwin/2; % time of the beginning of the spectrogram bin
end
toc
% % create matrix with zeroed signal
% templg = ones(2*n-1,length(S),'int8');
% templg(1,1:zwin) = 0; T = zeros(1,2*n-1);
% for i=2:n*2-1
%     templg(i,:) = circshift(templg(1,:),[0 zwin/2*(i-1)]);
%     T(i)=i*zwin/2;
% end
% gram = single(templg).*repmat(S,n*2-1,1);
% 
% % take powerspectrum of each zwin
% clear p;
% for i=1:size(gram,1)
%     p(i,:) = pwelch(gram(i,:),round(win),noverlap,nfft,1/dt);
% end
P = repmat(p0',size(p,1),1) - p;

%% plot
if bplot
figure(12);clf
ImageMatrix(T,F,temp')
axis xy;  colorbar('off'); 
ylim([0 200])
end
%