function [covxy time data] = plotxcov(x,y,range,dt,fnum,prange,varargin)
% function [covxy time data] = plotxcov(x,y,range,dt,fnum,prange)
% DESC:
% runs xcov(x,y,range,'coeff')
% shifts time axis to be centered at 0
% finds max correlation and time of max correlation
%
% INPUT
%   x, y:    vectors
%   range:   to run xcov over (samples)
%   dt:      sample time (s)
%   fnum:    number of figure to create
%   prange:  range to plot (ms) default [-100 100]
% runs xcov(x,y,range,'coeff')
% OUTPUT
%   covxy:    xcov
%   time:   lagtimes
%   data:      maxcor and maxtime (time of maxcor)
%
% Interpretation of Lag: if x PRECEDES y, the correlation will peak at a negative lag. 
% BA113006

%% default
scolor = '-k';
DIR = struct([]);
if nargin>=7
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

range = round(range);
bauto =0;
if isequal(size(x),size(y))
if sum(x-y)==0
    bauto = 1;
end
end
if ~exist('prange')
    prange = [-100 100];
end
d = find(size(x) ==min(size(x))); %get the right dimension
if d ~=1
    x = x';
    y = y';
end
covxy = multi_xcov(x,y,range)';
if size(x,1) >1
    covxy = nanmean(covxy,2);
end

if nargin < 5 || isempty(fnum)
    fnum = figure;
else
    figure(fnum);
end

[fnum time] = plotcorrelation(covxy,dt,fnum,prange,'color',scolor);

figure(fnum);
if ~bauto
    temp = covxy((size(covxy,1)-1))/abs(covxy((size(covxy,1)-1)));
    data.maxcor = max(abs(covxy))*temp;
    data.maxtime =  time(find(abs(covxy)==(max(abs(covxy)))));
    title(['CrossCov' ' peak: ' num2str(data.maxcor,'%1.2f') ' lag: ' num2str(data.maxtime,'%1.2f') 'ms']);
else
    title('AutoCorrelation')
end
plotset(1);