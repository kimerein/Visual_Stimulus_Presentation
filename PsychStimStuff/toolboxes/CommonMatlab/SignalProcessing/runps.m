function [p f param] = runps(x,dt,param)

if ~exist('dt')
    dt =1;
end
Res = 0.2;% (Hz)
win =[];
nwin = NaN;
noverlap = [];
nfft = 2^nextpow2((1/dt/2)/Res);
bdpss = 0;
if  nargin==3
    if ~isempty(param)
%         if isfield(param, 'novernfftlap')
        if isfield(param, 'nfft')
            nfft = param.nfft;
        else
            nfft = [];
        end
        if isfield(param, 'noverlap')
            noverlap =param.noverlap;
        end
        %     if isfield(param, 'nfftp')
        %         nfftp = param.nfftp;
        %     end
        if isfield(param, 'bdpss')
            bdpss = param.bdpss;
        end
        if isfield(param, 'nwin')
            nwin = param.nwin;
        end
        if isfield(param, 'win')
            win = param.win;
        end
    end
end
if isempty(win)
    if bdpss
        [win,v] = dpss(nwin,.5);
    end
end

if min(size(x))>1;  display('data must contain 1 trial per ROW');end
for i=1:size(x,1)
    [p(:,i) f] = pwelch(x(i,:),win,noverlap,nfft,1/dt);
end
if min(size(x))>1; p = mean(p,2); end;

param.win= win; param.nwin= nwin; param.bdpss = bdpss;
% param.nfftp= nfftp ;
param.noverlap= noverlap;
param.nfft       = nfft ;

