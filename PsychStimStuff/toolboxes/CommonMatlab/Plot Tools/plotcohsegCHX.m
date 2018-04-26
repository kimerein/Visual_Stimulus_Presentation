function [CxyCohere phi F S1 S2 confC phierr Cerr]= plotcohsegCHX(x,y,dt,win,varargin)
% function [CxyCohere phi F S1 S2 confC phierr Cerr]= plotcohsegCHX(x,y,dt,win,varargin)
% DESC: plots ouptput of coherencysegc.m (CHX)
% Note units have to be consistent. See chronux.m for more information.
%       x,y - data matrix each ROW new trial -- required
%       dt - sampling interval (default to 1)
%       win   (length of segments) - optional (default to entire length of
%       x)
%   varargin:
%          'fid':    number of figure to create
%          'fid_phase'
%          'color': color of line to plot
%           'bphase': boolean (1 plot phase) default 0
%           'bpw':   boolean (1 plot PS of x and y) default 0
%           'bplotlog:' coherence log plot
%       'params': structure with fields tapers, pad, fpass, err
%           tapers (precalculated tapers from dpss, or in the form [NW K] e.g [3 5]) -- optional. If not
%                                                 specified, use [NW K]=[3 5]
%	        pad		    (padding factor for the FFT) - optional. Defaults to 0.
%			      	 e.g. For N = 500, if PAD = 0, we pad the FFT
%			      	 to 512 points; if PAD = 2, we pad the FFT
%			      	 to 2048 points, etc.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional.
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%
% OUTPUT
%   CxyCohere:    real value
%   phi:   phase angle
%   F:     frequency
%
% BA061408
%% defaults
scolor = 'k';
bphase = 0;
bpw = 0;
cohparam.Fs = 1/dt;
cohparam.fpass = [0 200]; %kHz
cohparam.err = [0 0.05];
cohparam.tapers = [9 5];
bplotlog = 0;
fid = [];

if nargin <4 || isempty(win);    win = size(x,2); end
%% varargin
DIR = struct([]);
if nargin>=5
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
            case 'fid'
                fid = DIR(i).val;
            case 'bphase'
                bphase = DIR(i).val;
            case 'bpw'
                bpw = DIR(i).val;
            case 'cohparam'
                cohparam = DIR(i).val;
            case 'bplotlog'
                bplotlog = DIR(i).val;
            otherwise
        end
    end
end
% defaults for cohparam
cohparam.Fs = 1/dt;
if ~isfield(cohparam,'fpass'); cohparam.fpass = [0 200];  end
if ~isfield(cohparam,'err'); cohparam.err = [0 0.05]; end
if ~isfield(cohparam,'tapers'); cohparam.tapers = [3 5]; end
% if ~isfield(cohparam,'trialave'); cohparam.trialave = 1; end


if min(size(x))==1 & ~isrowvector(x)
    x = x';
end
if min(size(y))==1 & ~isrowvector(y)
    y = y';
end

Cerr = [];
for i = 1:size(x,1)
    if cohparam.err(1)==2
        [CxyCohere(:,i),phi(:,i),S12(:,i),S1(:,i),S2(:,i),F,confC,phierr,Cerr(:,:,i)]=coherencysegc(x(i,:)',y(i,:)',win,cohparam);
    elseif cohparam.err(1)
        [CxyCohere(:,i),phi(:,i),S12(:,i),S1(:,i),S2(:,i),F,confC,phierr]=coherencysegc(x(i,:)',y(i,:)',win,cohparam);
    else
        [CxyCohere(:,i),phi(:,i),S12(:,i),S1(:,i),S2(:,i),F]=coherencysegc(x(i,:)',y(i,:)',win,cohparam);
        Cerr = [];confC=[];phierr=[];
    end
end

if size(x,2) >1
    CxyCohere = mean(CxyCohere,2);
    Cerr = mean(Cerr,3); %% THIS may NOT be A GOOD IDEA!
    phi = mean(phi,2);
    S1 = mean(S1,2);
    S2 = mean(S2,2);
    %         S12 = mean(S12,2);
end
% if bphase | bpw
%
%     % Convert one-sided PSDs to Power Spectrums
%     Sxx = [Sxx(1); (Sxx(2:end-1)/2); Sxx(end)].*fs;
%     Syy = [Syy(1); (Syy(2:end-1)/2); Syy(end)].*fs;
%     % Cxy = (abs(Sxy).^2)./(Sxx.*Syy); %% to manually get amplitude
%     phase = unwrap(angle(Sxy)*180/pi);
% end
%%
if bplotlog; splt = 'l'; else splt = 'n'; end


if ~isempty(fid); figure(fid); else fid = gcf; end
% set(gca,'Position',[0 500 1000 1000],'Color','none','XGrid', 'on');
% subplot(1,2,1)
ax_h = plot_vector(CxyCohere',F,splt,Cerr,scolor);
stemp = [];
if cohparam.err(1); stemp = sprintf('  p > %1.1g',confC); end


% title(['Coherence' stemp]);
xlabel('freq (Hz)');
hold on;
plotset(1);
axis tight;
%
if bphase
    %     figure(fnum+1);
    % set(gca,'Position',[0 1500 1000 1000],'Color','none','XGrid', 'on');
    a = insetplot(.3,ax_h,1);
    phierr = []; % don't plot phierr
    plot_vector(phi,F,splt,phierr,scolor); % supposedly 95% confidence see coherencysegc.m

    ylabel('phase (degrees)')
    xlabel('')
    hold on;
    plotset(5);
    axis tight;
end

if bpw
    %     figure(fid+1);clf;
    % set(gca,'Position',[0 1500 1000 1000],'Color','none','XGrid', 'on');
    % subplot(1,2,2)
    insetplot(.3,ax_h,0);
    plot(F,S1/max(S1));hold all;
    plot(F,S2/max(S2));
    ylabel('Normalized Power spectrum abs^2/Hz')
    %     xlabel('freq (Hz)')
    xlim([0 100])
    hold on;
    plotset(1);
    axis tight;
end