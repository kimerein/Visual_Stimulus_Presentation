function playSTA(STA, STAparams, Nspikes, pausetime)
% function playSTA(STA, STAparams,pausetime)
%

if nargin < 3; pausetime = 0.2; end

tbin = STAparams.tbin;
WOI = STAparams.WOI;
figure(1);
for i=1: size(STA,3)
    %     subplot(2,1,1)
    imagesc(STA(:,:,i));
    %     subplot(2,1,2)
    %     imagesc(wSTA(:,:,i));
    title(['spk:' num2str(Nspikes) ' ' num2str(i) ' ' num2str((i*tbin-WOI(1))*1000) 'ms']);
    pause(pausetime)
end