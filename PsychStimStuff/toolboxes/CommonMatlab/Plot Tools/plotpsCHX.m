function [S f Serr] = plotpsCHX(x,dt,fid,prange,param,pparam,bpt)
% function [p1 f] = plotpsCHX(x,dt,fid,prange,param,pparam)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data (in form samples x channels/trials) -- required
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%       -optional
%           tapers (precalculated tapers from dpss, or in the form [NW K] e.g [3 5]) -- optional. If not 
%                                                 specified, use [NW K]=[3 5]
%	        pad		    (padding factor for the FFT) - optional. Defaults to 0.  
%			      	 e.g. For N = 500, if PAD = 0, we pad the FFT 
%			      	 to 512 points; if PAD = 2, we pad the FFT
%			      	 to 2048 points, etc.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%           trialave (average over trials/channels when 1, don't average when 0) - optional. Default 0
% Output:
%       S       (spectrum in form frequency x channels/trials if trialave=0; in the form frequency if trialave=1)
%       f       (frequencies)
%       Serr    (error bars) only for err(1)>=1
% BA 060208 functino to plot from Chronux
% pparam.logplot = default 0
if nargin <3 || isempty(fid)
 fid=    figure;
else
    figure(fid);
end

if nargin <4 || isempty(prange)
    prange = [0 100];
end
if nargin <7; bpt=0; end

% defaults
iparam.Fs = 1/dt;
iparam.trialave = 1; % average sweeps

if nargin >=5;
    if  isfield(param, 'trialave');        iparam.trialave = param.trialave;    end
    if  isfield(param, 'err');        iparam.err = param.err;    end
    if  isfield(param, 'Fs');        iparam.Fs = param.Fs;    end
    if  isfield(param, 'fpass');        iparam.fpass = param.fpass;  prange = iparam.fpass;  end
end

bploterr = 0; % plot error bars
if isfield(iparam,'err')
    if iparam.err(1)>0;   bploterr=1; end
end

st = 'n'; c= [0 0 0];
if  exist('pparam')
    if isfield(pparam, 'logplot')
        if pparam.logplot; st = 'l'; end
    end
    if isfield(pparam, 'color')
        c= pparam.color;  
    end
end

if bploterr
    if ~bpt
        printf('Continuous signal PS')
        [S,f,Serr]=mtspectrumc(x,iparam);
    else
        [S,f,Serr]=mtspectrumpb(x,iparam);
    end
    plot_vector(S,f,st,Serr,c)
else
    if ~bpt
                printf('Continuous signal PS')
        [S,f]=mtspectrumc(x,iparam);
    else
        [S,f]=mtspectrumpb(x,iparam);
    end
    
    plot_vector(S,f,st,[],c)
    Serr = [];
end
plotset(1);
xlim(prange)
% title(['PS  win:' num2str(param.nwin*dt,'%1.2f') 's nfft:' num2str(param.nfft*dt,'%1.2f') 's']);
xlabel('Frequency (Hz)')
ylabel('(units^2/Hz)')
