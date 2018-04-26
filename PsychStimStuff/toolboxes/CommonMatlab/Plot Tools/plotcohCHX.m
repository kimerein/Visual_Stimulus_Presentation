function [CxyCohere phi F S1 S2 confC phierr Cerr]=  plotcohCHX(x,y,dt,varargin)
% function [CxyCohere phi F S1 S2 confC phierr Cerr]= plotcohCHX(x,y,dt,varargin)
% DESC: plots ouptput of coherencyc.m (CHX)
% Note units have to be consistent. See chronux.m for more information.
%       x,y - data matrix each ROW new trial -- required
%       dt - sampling interval (default to 1)
%       x)
%   varargin: 
%          'fid':    number of figure to create
%          'color': color of line to plot
%           'bphase': boolean (1 plot phase) default 0
%           'bpw':   boolean (1 plot PS of x and y) default 0
%           'bplotlog:' coherence log plot
%           'segwin' window size to segment data before doing coherence (in units of dt)
%                    used in coherencysegc, Note: if segwin is not defined
%                    then coherencyc (i.e. no segmentation of time series)
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
scolor = '-k';
bphase = 0;
bpw = 0;
bplotlog = 0;
segwin = []; % 
fid = [];

if nargin <4 ;    win = size(x,2); end
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
            case 'segwin'
                segwin = DIR(i).val;
            otherwise
        end
    end
end

% defaults for cohparam
cohparam.Fs = 1/dt;
if ~isfield(cohparam,'fpass'); cohparam.fpass = [0 200];  end
if ~isfield(cohparam,'err'); cohparam.err = [0 0.05]; end
if ~isfield(cohparam,'tapers'); cohparam.tapers =[3 5]; end
if ~isfield(cohparam,'trialave'); cohparam.trialave = 1; end

if min(size(x))==1 & ~isrowvector(x)
    x = x';
end
if min(size(y))==1 & ~isrowvector(y)
    y = y';
end

if cohparam.err(1)
    if isempty(segwin);        [CxyCohere,phi,S12,S1,S2,F,confC,phistd,Cerr]=coherencyc(x',y',cohparam);
    else;        [CxyCohere,phi,S12,S1,S2,F,confC,phistd,Cerr]=coherencysegc(x',y',segwin,cohparam);    end
    % it seems easiest to control bandwidth with limiting the window size
    % using this function

else

    if isempty(segwin);        [CxyCohere,phi,S12,S1,S2,F]=coherencyc(x',y',cohparam);
    else [CxyCohere,phi,S12,S1,S2,F]=coherencysegc(x',y',segwin,cohparam);    end

    Cerr = [];        Cerr = [];confC=[];phistd=[];

end

%%
if bplotlog; splt = 'l'; else splt = 'n'; end


if ~isempty(fid); figure(fid); end
% set(gca,'Position',[0 500 1000 1000],'Color','none','XGrid', 'on');
if bphase;subplot(1,2,1);end
ax_h = plot_vector(CxyCohere,F,splt,Cerr,scolor); % convert to Hz

stemp = [];
if cohparam.err(1); stemp = sprintf('  p > %1.1g',confC); end


title(['Coherence' stemp]);
xlabel('freq (Hz)');
hold all;
plotset(1);
%
if bphase
%     figure(fnum+1);
    % set(gca,'Position',[0 1500 1000 1000],'Color','none','XGrid', 'on');
    subplot(1,2,2)
%     insetplot(.3,ax_h,1);
%     phistd = []; % don't plot phistd
 plot_vector(phi,F,splt,phistd,scolor);
    ylabel('phase (degrees)')
    xlabel('freq (Hz)')
    hold all;
    plotset(5);
end

if bpw
    figure(fnum+1);clf;
    % set(gca,'Position',[0 1500 1000 1000],'Color','none','XGrid', 'on');
    % subplot(1,2,2)
       insetplot(.3,ax_h,0);

    plot(F,S1/max(S1));hold all;
    plot(F,S2/max(S2));
    xlim(prange) ;
    ylabel('Normalized Power spectrum abs^2/Hz')
    xlabel('freq (Hz)')
    hold on;
    plotset(1);
end