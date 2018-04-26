function h = plotAnn(stmp,fid,pos,sTag)
% function h = plotAnn(stmp,fid,pos,sTag)
% Annotate plot: places text in the left margin of the figure

if nargin >1
    if  ~isempty(fid)
        set(0,'CurrentFigure',fid);
    else
        
    end
else
    fid = gcf;
end
if isempty(stmp)
    stmp = '';
end
ax = gca;
if nargin<4; sTag = ''; end
% Ccolor = get(gcf,'Color');a
Ccolor = 'None';
if nargin<3 || isempty(pos)
    pos = 3;
end
hold on
% stmp = 'test';
switch(pos)
    case 1
        axes('position',[0 1-.05  1 .04]);
        h = text(1-0.02,.5,stmp,'rotation',0,'Interpreter','none','HorizontalAlignment','right','BackgroundColor',Ccolor,'Tag',sTag);
    case 2
        axes('position',[1-.02 0.02  0.02 1]);
        h = text(0.02,0,stmp,'rotation',90,'Interpreter','none','HorizontalAlignment','left','BackgroundColor',Ccolor,'Tag',sTag);
    case 3
        axes('position',[0.01 0.02  1-0.01 0.02]);
        h = text(1-.01,0.02,stmp,'rotation',0,'Interpreter','none','HorizontalAlignment','right','BackgroundColor',Ccolor,'Tag',sTag);
    case 4
        axes('position',[0 .05  .04 1]);
        h = text(.5,0,stmp,'rotation',90,'Interpreter','none','BackgroundColor',Ccolor,'Tag',sTag);
end
axis off;
set(fid,'CurrentAxes',ax);