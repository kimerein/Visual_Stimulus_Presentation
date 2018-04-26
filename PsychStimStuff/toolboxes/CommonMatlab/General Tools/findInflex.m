function [indInfx bpeak] = findInflex(d);
% function [indInfx bpeak] = findInfx(d);
% Find inflection points (i.e. find peaks and troughs) in d
% classify as peak (or trough)
%
% INPUT: d: data vector
%  
% OUTPUT 
%        indInf - indices of inflection points in d
%   NOTE: since data is discrete and there may not actually be an zero value in
%   the derivative. function returns the index of the last value before
%   zero is crossed.
%        bpeak - 1 if peak 0 if trough
%   NOTE: technically only looks at the slope of the derivative before the
%   zero crossing (not both before and after) so can't tell if there were a
%   real inflection point.
% BA012907
%
% Useful: for analyzing gamma oscillations

% find zero crossings
dd = single(diff(d));
d2d = single(diff(dd));
% if the product of two sequential values is <0 there was a zero crossing
% (peak or trough)
tmp = circshift(dd,1);
tmp = dd(2:end).* tmp(2:end);
indInfx = find(tmp <0) + 2; % index of peaks and trough (add 2 to correct for 2 shifts in index compared to original data)

bpeak = zeros(size(indInfx),'single');
bpeak(find(d2d(indInfx-2)<0)) = 1;% classify as peak or trough (peak if 2nd Derivitive is < 0)