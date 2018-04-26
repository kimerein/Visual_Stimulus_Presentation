function [th ind_xth sw] = threshdata(data,dt,th,bplot,minSampleAboveThres)
% function [th ind_xth sw] = threshdata(data,dt,th,bplot)
%
% OUTPUT: th threshold
%          ind_xth indices of threshold crossing in data
%
% !!NOTE: threshold looks for events below th
% BA
if ~exist('minSampleAboveThres'); minSampleAboveThres = []; end
if nargin <3 | isempty(th)
    th = UIgetThresh(data,dt);
end
if nargin <4 ; bplot = 1; end


if nargout >2
    [ind_xth sw] = thresh(data,th,minSampleAboveThres);
else
    ind_xth = thresh(data,th);
end
if bplot
    plotev(data,ind_xth,dt,th,200);
end