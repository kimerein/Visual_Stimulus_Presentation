function formPlotwv(loc,lb)
% function formPlotwv(loc,lb)
% formats current axis for a waveform
% loc specifies the location of scale bars
%   1  2
%   3  4
% lb is a cell that specifies units for the scalebar 
%  e.g. lb{1} = 'mV' (y-label)
%  lb{2} = 'ms' (x-label)
%
% BA 101407
if nargin < 2
    lb{1} = '';
    lb{2} = '';
elseif size(lb)<2
    lb{2} = '';
end
set(gca,'Color','none')
a1 = gca;
axis off

%%
w = 0.01; % width of bars
%
a = axis();
xr = range([a(1),a(2)]);
x = xr*.1; % bar should be tenth of full range
i=0;
while x>10 %round to the nearest most signifigant digit
    i = i+1;
    x = x/10;
end
x = round(x)*10^i;
%
yr = range([a(3),a(4)]);
y = yr*.1; % bar should be tenth of full range
i=0;
while y>10 %round to the nearest most signifigant digit
    i = i+1;
    y = y/10;
end
y = round(y)*10^i;
%%
p = get(gca,'Position');
a2 = axes('Position',p);
hold on
axis off

% X bar%\
clear st;
switch(loc)
    case 1
        rectangle('Position',[0,1-w,x/xr,w],'FaceColor','k')
        rectangle('Position',[0,1-y/yr,w,y/yr],'FaceColor','k')
        textpos(1) = w*3; textpos(2) = 1-w*6; textalgn = 'left';
    case 2
        rectangle('Position',[1-x/xr,1-w,x/xr,w],'FaceColor','k')
        rectangle('Position',[1-w,1-y/yr,w,y/yr],'FaceColor','k')
       textpos(1) = 1-w*3; textpos(2) = 1-w*6; textalgn = 'right';
 
      case 3
        rectangle('Position',[0,0,x/xr,w],'FaceColor','k')
        rectangle('Position',[0,0,w,y/yr],'FaceColor','k')
      textpos(1) = w*3; textpos(2) = w*6; textalgn = 'left';
 
     case 4
        rectangle('Position',[1-x/xr,0,x/xr,w],'FaceColor','k')
        rectangle('Position',[1-w,0,w,y/yr],'FaceColor','k')
      textpos(1) = 1-w*3; textpos(2) = w*6; textalgn = 'right';
 
end
        st(1) = {[num2str(y) ' ' lb{1}]};
        st(2) = {[num2str(x) ' ' lb{2}]};

        text(textpos(1),textpos(2),st,'HorizontalAlignment',textalgn)

axis([0 1 0 1])
set(gca,'Color','none')
axes(a1);

