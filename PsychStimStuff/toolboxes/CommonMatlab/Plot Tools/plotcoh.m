function [CxyCohere F1 phase F Sxx Syy] = plotcoh(x,y,dt,fnum,prange,bphase,param,bpw,varargin)
% function [CxyCohere F1 phase F Sxx Syy]= plotcoh(x,y,dt,fnum,prange,bphase,param,bpw)
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

scolor = '-k';
DIR = struct([]);
if nargin>=9
    for i=1:length(varargin)
        if mod(i,2)~=0
        DIR(floor(i/2)+1).param = varargin{i};
        else 
        DIR(floor(i/2)).val = varargin{i};
        end
    end
end
for i=1:length(DIR)
    if ~isempty(DIR(i).param)&~isempty(DIR(i).val)
        switch DIR(i).param
            case 'color'
                scolor = DIR(i).val;
            otherwise
        end
    end
end

if nargin < 4 || isempty(fnum)
    fnum = 3000;
end
if nargin <5 || isempty(prange)
    prange = [0 100];
end
if nargin <6 || isempty(bphase)
    bphase = 1;
end
if nargin <8
    bpw= 1;
end


phase = []; F = [];

nfft = 8*2048;
nfftp = 4*2048;
noverlap = 2^8;
win = [];hamming(2^12);
bdpss = 0;
noverlapp = [];
%% parameters
if  exist('param')
    if isfield(param, 'novernfftlap')
        nfft = param.nfft;
    else
        nfft = [];
    end
    if isfield(param, 'noverlap')
        noverlap =param.noverlap;
    else
        noverlap = [];
    end
    if isfield(param, 'nfftp')
        nfftp = param.nfftp;
    else
        nfftp = [];
    end
    if isfield(param, 'bdpss')
        bdpss = param.bdpss;
    end
    if isfield(param, 'nwin')
        nwin = param.nwin;
    else
        nwin = [];
    end
    if isfield(param, 'win')
        win = param.win;
    end
end
% if 1/dt > length(x) |1/dt > length(y)
%     nwin = hamming(2^nextpow2(min(length(x),length(y))));
% else
%     nwin = hamming(2^nextpow2(1/dt)); % choose resolution of 1Hz (unless x or y are too short)
% end
if isempty(win)
    if bdpss
        [win,v] = dpss(nwin,.5);
    else
        win = hamming(min([2^12,length(x),length(y)])) ;% default
    end
end

nwinp = hamming(min([2^10,length(x),length(y)])); % default
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
    [CxyCohere(:,i),F1] = mscohere(detrend(x(:,i)),detrend(y(:,i)),win,noverlap,nfft,fs);

    %% Get Phase
    % Compute the PSDs using Welch's method, which is similar to the
    % MSCOHERE function
    if bphase |bpw

        [Sxx(:,i) ,F] = pwelch(x(:,i),nwinp,noverlapp,nfftp,fs);
        [Syy(:,i) ,F] = pwelch(y(:,i),nwinp,noverlapp,nfftp,fs);
        % Call CPSD function with similar arguments as the MSCOHERE
        % function from above.
        [Sxy(:,i) ,F] = cpsd(x(:,i),y(:,i),nwinp,noverlapp,nfftp,fs);
    end
end

if size(x,2) >1
    CxyCohere = mean(CxyCohere,2);
    if bphase | bpw

        Sxx = mean(Sxx,2);
        Syy = mean(Syy,2);
        Sxy = mean(Sxy,2);
    end
end
if bphase | bpw

    % Convert one-sided PSDs to Power Spectrums
    Sxx = [Sxx(1); (Sxx(2:end-1)/2); Sxx(end)].*fs;
    Syy = [Syy(1); (Syy(2:end-1)/2); Syy(end)].*fs;
    % Cxy = (abs(Sxy).^2)./(Sxx.*Syy); %% to manually get amplitude
    phase = unwrap(angle(Sxy)*180/pi);
end
%% Plot
figure(fnum);
% set(gca,'Position',[0 500 1000 1000],'Color','none','XGrid', 'on');
% subplot(1,2,1)
plot(F1, [CxyCohere],scolor);
title('Coherence');
%     plot(F1, Cxy); %% can alternatively use manually computed value
%     (similar but not identical)
xlim(prange);
xlabel('freq (Hz)');
hold on;
plotset(1);
%
if bphase
%     figure(fnum+1);
    % set(gca,'Position',[0 1500 1000 1000],'Color','none','XGrid', 'on');
    % subplot(1,2,2)
    insetplot(.3,[],1);
    plot(F,phase,scolor);
    xlim(prange) ;
    ylabel('phase (degrees)')
    xlabel('freq (Hz)')
    hold on;
    plotset(5);
end

if bpw
    figure(fnum+2);clf;
    % set(gca,'Position',[0 1500 1000 1000],'Color','none','XGrid', 'on');
    % subplot(1,2,2)
    plot(F,Sxx/max(Sxx));hold all;
    plot(F,Syy/max(Syy));
    xlim(prange) ;
    ylabel('Normalized Power spectrum abs^2/Hz')
    xlabel('freq (Hz)')
    hold on;
    plotset(1);
end