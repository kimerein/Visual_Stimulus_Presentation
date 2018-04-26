function xcorrmk(X,Y,allspikeTimes,ID,Ts,varargin)
% function  xcorrmk(X,Y,allspikeTimes,ID,Ts,varargin)
%
% DESC: crosscorrelated spikes with labels in X with those with labels in Y.
% used to cross-correlated isolated units or klusters
%
% INPUT:
%       X - vector of labels (e.g. [1 2 4 5..]
%       Y - vector of labels
%       allspikeTimes - time of spikes in experiment (samples)
%       ID - label for each of the spikes in allspikeTimes, size(ID)==size(allspikeTimes)
%       Ts - sampling interval (sec)
% optional:
%        fid - figure id
%        range - in ms of cross correlation [500]
%          specify as 'range',range
%        bsz - in ms [1]
%          specify as 'bsz',bsz
%        color - 
%
% OUTPUT: length(X) figures each with length(Y) subplots
%         can also save to disk
%
% NOTE: xcorr is done manually and range is hardcoded
%                                    binsize is hardcoded
% BA 012307
%defaults
% srebin = 1;
% range =500;
srebin = 0.5;
range =500;
fid = 700;

DIR = struct([]);
if nargin>=6
    for i=1:length(varargin)
        if mod(i,2)~=0
        DIR(floor(i/2)+1).param = varargin{i};
        else 
        DIR(floor(i/2)).val = varargin{i};
        end
    end
end

scolor = '-b';
for i=1:length(DIR)
    if ~isempty(DIR(i).param)&~isempty(DIR(i).val)

        switch DIR(i).param
            case 'range'
                range = DIR(i).val;
            case 'bsz'
                srebin = DIR(i).val;
            case 'fid'
                fid = DIR(i).val;
            case 'color'
                scolor = DIR(i).val;
            otherwise
        end
    end
end


range =range/srebin;
binsize = srebin*1e-3;
binperdata = Ts/1000;
figheader = 'KlustersXcorr';

NX = size(X,2);
NY = size(Y,2);
NUnits = NY;
NS = min(5,NUnits);
COL = ceil(NUnits/NS);
NN = NX;
for k = 1:NN % Number of figures (1 for each cell)
    figure(k-1+fid);
%       set(gcf,'Position',[ 2635        -141         643        1120])

    hold on
      
    unitX = X(k);
    unit_spiketimesX = allspikeTimes(find(ID == unitX));
    
            binned_spikesX = rebin(bin2(binned_spikesX,binperdata)',srebin);

%     binned_spikesX = rasterize(unit_spiketimesX,srebin);
   %% rebin normalization so that area is remains correct 
    for j = k:NY %% Number of subplots per figure (1 for each cell)
    unitY =Y(j);
    unit_spiketimesY = allspikeTimes(find(ID == unitY));
    
    if NX>1 && NY>1
        subplot(NS,COL,j)
    end
       % Bin Data
        binned_spikesY = rebin(bin2(unit_spiketimesY,binperdata)',srebin);
%         binned_spikesY =rasterize(unit_spiketimesY,srebin);
        %         plot((1:size(binned_spikes,2))*(binsize),binned_spikes(1,:))  % X axis time in Se

        cross_corr = xcorrmanual(binned_spikesX,binned_spikesY,range);
        xtime = ([1:size(cross_corr,2)]-(range+1))*binsize*1000;
        if ~isstr(scolor)
            plot(xtime,cross_corr,'color',scolor);% X axis time in Sec
        else
            plot(xtime,cross_corr,scolor);% X axis time in Sec
        end
%         plot(cross_corr,'-b');% X axis time in Sec
        axis tight
        xlabel('ms');
%         temp = sprintf('%s%d-%d (%d,%d) Bin:%1.2gms',figheader,k,j,unit_Nspikes(k),unit_Nspikes(j),binsize*1000);
        temp = sprintf('%s%d-%d Spikes:%d,%d*, Bin:%1.2gms',figheader,unitX,unitY,length(unit_spiketimesX),length(unit_spiketimesY),binsize*1000);
        title(temp);
    end
%     %% save data image to disk
%     dirtemp = 'REPORT';
%     figdesc = 'ManualXCorr';
%     spath = [writedirheader dirtemp '\'];
%     if (~isdir(spath))
%         mkdir(spath);
%     end
%     temp = sprintf('%s%s_Unit_%d.png',spath,figdesc,unitX);
% %     print('-dtiffn',temp);
%     print('-dpng',temp);
end