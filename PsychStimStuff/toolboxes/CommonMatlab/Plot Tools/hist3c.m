function [n C] =hist3c(x,y, nbins,fid,bplotlog)
% function [n C] =hist3c(x,y, nbins,fid,bplotlog)
%
% takes same input as hist3.
% creates a 3d histogram of x and y using bins(1) for x and nbins(2) (default 10) for y
% fid - figure handle (optional)
%
if nargin <5
bplotlog = 0;
end
if nargin <3 || isempty(nbins)
    nbins = [10 10];
end

if isrowvector(x)
    x = x';
end
if isrowvector(y)
    y = y';
end

[n C] = hist3([y x],nbins); % flip order of x and y so it comes out like traditional plot

if nargin <4 || isempty(fid)
    figure;
else
    figure(fid);
end

if bplotlog
    n = log2(n);
end
imagesc(C{2},C{1},-n)
axis xy

axis tight
colorbar
colormap(gca,'gray')