function plotscalebar(loc,lb)
% function plotscalebar(loc,lb)
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
elseif sum(size(lb))<2
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
if x<1
    while x<1
        i = i+1;
        x = x*10;
    end
    x = round(x)*10^-i;
elseif x>10

    while x>10 %round to the nearest most signifigant digit
        i = i+1;
        x = x/10;
    end
    x = round(x)*10^i;
else
    x = round(x);
end
%
yr = range([a(3),a(4)]);
y = yr*.1; % bar should be tenth of full range
i=0;
if y<1
    while y<1
        i = i+1;
        y = y*10;
    end
    y = round(y)*10^-i;
elseif y>10

    while y>10 %round to the nearest most signifigant digit
        i = i+1;
        y = y/10;
    end
    y = round(y)*10^i;
else
    y = round(y);
end
%%
p = get(gca,'Position');
a2 = axes('Position',p);
axis([0 1 0 1])

hold on
axis off
LINEWIDTH = 3;
% X bar%\
clear st;
switch(loc)
    case 1
        %         rectangle('Position',[0,1-w,x/xr,w],'FaceColor','k')
        %         rectangle('Position',[0,1-y/yr,w,y/yr],'FaceColor','k')
        %                 textpos(1) = w*3; textpos(2) = 1-w*11; textalgn =
        %                 'right';
        line([0 x/xr],[1 1],'Color','k','LineWidth',LINEWIDTH)
        line([0,0],[1-y/yr,1],'Color','k','LineWidth',LINEWIDTH)
        textpos(1) = w*3; textpos(2) = 1 ;textalgn = 'left';
        st(2) = {[num2str(y) ' ' lb{1}]};
        st(1) = {[num2str(x) ' ' lb{2}]};

    case 2
        %         rectangle('Position',[1-x/xr,1-w,x/xr,w],'FaceColor','k')
        %         rectangle('Position',[1-w,1-y/yr,w,y/yr],'FaceColor','k')
        line([1-x/xr 1],[1 1],'Color','k','LineWidth',LINEWIDTH)
        line([1,1],[1-y/yr,1],'Color','k','LineWidth',LINEWIDTH)
        textpos(1) = 1-w*3; textpos(2) = 1; textalgn = 'right';
        st(2) = {[num2str(y) ' ' lb{1}]};
        st(1) = {[num2str(x) ' ' lb{2}]};
    case 3
        %         rectangle('Position',[0,0,x/xr,w],'FaceColor','k')
        %         rectangle('Position',[0,0,w,y/yr],'FaceColor','k')
        %                 textpos(1) = w*3; textpos(2) = w*11; textalgn = 'left';

        line([0 x/xr],[0 0],'Color','k','LineWidth',LINEWIDTH)
        line([0,0],[0,y/yr],'Color','k','LineWidth',LINEWIDTH)
        textpos(1) = w*3; textpos(2) = 0; textalgn = 'left';
        st(1) = {[num2str(y) ' ' lb{1}]};
        st(2) = {[num2str(x) ' ' lb{2}]};

    case 4
        %         rectangle('Position',[1-x/xr,0,x/xr,w],'FaceColor','k')
        %         rectangle('Position',[1-w,0,w,y/yr],'FaceColor','k')
        %                 textpos(1) = 1-w*3; textpos(2) = w*11; textalgn =
        %                 'right';
        line([1-x/xr 1],[0 0],'Color','k','LineWidth',LINEWIDTH)
        line([1,1],[0,y/yr],'Color','k','LineWidth',LINEWIDTH)
        textpos(1) = 1-w*3; textpos(2) = 0; textalgn = 'right';
        st(1) = {[num2str(y) ' ' lb{1}]};
        st(2) = {[num2str(x) ' ' lb{2}]};

end
% st(1) = {[num2str(y) ' ' lb{1}]};% for rectangle
% st(2) = {[num2str(x) ' ' lb{2}]};

text(textpos(1),textpos(2),st,'HorizontalAlignment',textalgn)

axis([0 1 0 1])
set(gca,'Color','none')
axes(a1);

