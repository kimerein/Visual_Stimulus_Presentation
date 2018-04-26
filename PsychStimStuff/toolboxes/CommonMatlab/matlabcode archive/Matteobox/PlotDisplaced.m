function Displacement = PlotDisplaced(tt,zz,v,mode)
% PlotDisplaced plots traces vertically displaced
% 
% PlotDisplaced(tt,zz) with zz nTraces X nSamples
%
% PlotDisplaced(tt,zz,v) lets you specify how many standard deviations v
% separate the traces (DEFAULT: 5)
%
% PlotDisplaced(tt,zz,v,'absolute') interprets v as an absolute displacement (not
% the number of standard deviations)
%
% Displacement = PlotDisplaced(...) reports the absolute displacement
%
% 2007 Matteo Carandini
% 2007-06 SK removed limitation on number of traces by looping through 'colors' string 
% 2007-06 MC added an indication of amplitude
% 2007-06 MC added output of displacement

if nargin<4
    mode = 'std';
end

if nargin<3
    v = 5;
end

[nc,nt] = size(zz);

%colors = 'rbgkrbgkrbgkrbgkrbgkrbgk';
colors = 'rbgk';

switch mode
    case 'std'
        Displacement = v*nanstd(zz(:));
    otherwise
        Displacement = v;
end

for ic = 1:nc
    plot( tt, zz(ic,:) - (ic-1)*Displacement, 'color', colors(mod(ic-1, length(colors))+1) ); hold on
end
axis tight
set(gca,'ytick',zz(1,1)+[-Displacement 0 ],'yticklabel',[0 str2num(num2str(Displacement,2))],'ycolor','k','box','off');

