function a = insetplot(frac,axishandle,corner)
% function a = insetplot(frac,axishandle,corner)
% creates plot in top right hand, with dimsions fraction, frac, of plot.
% INPUT:
%        frac: fraction of current axis to make inset.
% NOTE: plot dimensions are square, defined frac*min(width,height)
%        axishandle: (optional) handle of axis to use default current axis
%        corner: (optional)  0 topleft 1 topright,
% EXAMPLE: for an inset that is 20% the size of current plot use:
%          insetplot(0.2)
%This is written to place inset in top right corner but can easily be
%expanded
%
%BA050108
if (nargin >1 && ~isempty(axishandle)); axes(axishandle); end
if ( nargin < 2); corner=1; end
% aa = gca;
Xedgeoffset = 0.009;
Yedgeoffset = 0.009;

p = get(gca,'Position');
% plot dimensions are square, defined frac*min(width,height)
W = frac*min(p(3:4));
H = W;
switch (corner)
    case 1
        Xpos = (p(1) + p(3)) - Xedgeoffset - W;
        Ypos = (p(2) + p(4)) - Yedgeoffset - H;
    case 0
        Xpos = p(1) + Xedgeoffset*2;
        Ypos = (p(2) + p(4)) - Yedgeoffset - H;
end

axes('Position',[Xpos Ypos W H],'box','off');
a = gca;
% axes(aa); % return focus to last axis