function [fnum time] = plotcorrelation(x,dt,fnum,prange,varargin)
% function [fnum time] = plotcorrelation(x,dt,fnum,prange,varargin)
% function to help plot correlations data so that x axis is centered at
% zero lag
%
% OUTPUT fnum = figure
%        time = values of x axis (ms if dt is in secounds)
%
% BA051508
scolor = '-';
DIR = struct([]);
if nargin>=4
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

time = plotcorrelationHelper(size(x,1),dt);

if nargin < 3 || isempty(fnum)
    figure;
else
    figure(fnum);
end

plot(time,x,scolor,'LineWidth',2)
hold on
axis tight
if exist('prange'); xlim(prange); end
xlabel('time (ms)')
