function handles = scope(name, data, xData)
% function handles = scope(name, data, xData)
% display data channels in a scrollable window
% scope(chanData)
% scope(name, chanData)
% scope(chanData, timeData)
% scope(name, chanData, timeData)
%
% left click and drag on the axes to zoom in, right click to zoom out.
% clicking on the display will set or unset the cursors.
% the text box of the bottom left corner is the zoom factor
% scroll up and down to change sweeps % BA
%
% TO DO: ADD    Display sweep number
%               manually define range of sweeps,
%               speed up display of plots
% parse input
switch nargin
    case 1
        data = name;
        name = 'Scope';
        xData = 1:length(data);
    case 2
        if ischar(name)
            xData = 1:length(data);
        else
            xData = data;
            data = name;
            name = 'Scope';
        end
    case 3
        % do nothing
    otherwise
        disp('Improper input')
end

if size(data, 1) < size(data, 2)
    data = data';
end

if size(data) == 1
    error('data is size 1') % BA
    %     noData = 1;
    %     data = ones(1,data);
    % else
    %     noData = 0;
end

%initialize plotting window
set(0,'Units','normal');

figure('NumberTitle','off',...
    'Name', name,...
    'menu', 'none',...
    'Units','normal',...
    'Position',[0 .025 1 .95],...
    'Visible', 'on',...
    'UserData', -100,...
    'windowButtonDownFcn', @mouseDownScope,...
    'WindowButtonMotionFcn', @movePointers,...
    'windowButtonUpFcn', @mouseUpScope,...
    'WindowScrollWheelFcn',@mousescroll); % BA
clear control_info;

uicontrol('Style','slider','Units','normal','Position', [0 0 1 .015], 'value', 0, 'sliderStep', [1 1/0], 'callback', @traceScroll);
uicontrol('Style','edit','Units','normal','Position',[0 .018 .02 .02], 'string', '1', 'callback', @traceScroll);

sw = 1;
global SCOPExData SCOPEdata SCOPEsw SCOPEhandles

SCOPExData = xData; SCOPEdata = data; SCOPEsw = sw;
for index = 1:size(SCOPEdata, 2); % nchns
    SCOPEhandles(index) = axes('Position', [.05 .05 + (index - 1) / (size(SCOPEdata, 2) / .95)  .95 .95 / size(SCOPEdata, 2)], 'nextplot', 'add');
    startLine = line([1 1], [0 1], 'color', [0 0 0], 'userData', 0);
    stopLine = line([0 0], [0 1], 'color', [0 0 0]);
    lineText = text(1, 1, '0', 'color', [0 0 0], 'VerticalAlignment', 'bottom', 'linestyle', 'none');
    text(0, 1, '0', 'visible', 'off', 'VerticalAlignment', 'bottom', 'color', [0 0 0], 'linestyle', 'none')
    if index > 1
        set(gca, 'xticklabel', '', 'xtick', []);
    end
end
plotdatahelper();


set(gcf, 'userData', 0);

function mouseDownScope(varargin)
pointerLoc = get(gcf, 'CurrentPoint');
imageLoc = get(gca, 'Position');
if pointerLoc(1) > imageLoc(1) && pointerLoc(1) < imageLoc(1) + imageLoc(3) && pointerLoc(2) > .05
    switch get(gcf, 'SelectionType')
        case 'normal' %left mouse button clicked
            kids = get(gcf, 'children');
            firstKid = get(kids(1), 'children');
            if get(firstKid(length(firstKid)), 'userData') == 0
                for index = 1:length(kids) - 2
                    kidLines = get(kids(index), 'children');
                    set(kidLines(length(kidLines)), 'color', [1 0 0]);
                end
                set(firstKid(length(firstKid)), 'userData', 1);
            else
                for index = 1:length(kids) - 2
                    kidLines = get(kids(index), 'children');
                    set(kidLines(length(kidLines)), 'color', [0 0 0]);
                    set(kidLines(length(kidLines) - 1), 'xData', [min(get(kids(index), 'xlim')) min(get(kids(index), 'xlim'))]);
                    set(kidLines(length(kidLines) - 3), 'visible', 'off')
                end
                set(firstKid(length(firstKid)), 'userData', 0);
            end
        case 'extend' %middle mouse button clicked

        case 'alt' %right mouse button clicked

        case 'open' %double click

    end
elseif pointerLoc(1) <= imageLoc(1) && pointerLoc(2) > 0.05
    kids = get(gcf, 'children');
    whichAxis = fix((pointerLoc(2) - .05) / .95 * (length(kids) - 2)) + 1;
    switch get(gcf, 'SelectionType')
        case 'normal' %left mouse button clicked
            set(gcf, 'userData', [whichAxis pointerLoc(2)]);
            uicontrol('units', 'normal', 'Style', 'text', 'String', '', 'backgroundColor', [0 0 1], 'Position', [0, pointerLoc(2), .01, 0.0001]);
        case 'extend' %middle mouse button clicked

        case 'alt' %right mouse button clicked

        case 'open' %double click

    end
elseif pointerLoc(1) > imageLoc(1) && pointerLoc(1) < imageLoc(1) + imageLoc(3) && pointerLoc(2) <= .05
    %kids = get(gcf, 'children');
    switch get(gcf, 'SelectionType')
        case 'normal' %left mouse button clicked
            set(gcf, 'userData', [-1 pointerLoc(1)]);
            uicontrol('units', 'normal', 'Style', 'text', 'String', '', 'backgroundColor', [0 0 1], 'Position', [pointerLoc(1), .015, 0.0001, .01]);
        case 'extend' %middle mouse button clicked

        case 'alt' %right mouse button clicked

        case 'open' %double click

    end
end

function mousescroll(obj,scrolldata)
global SCOPEsw
SCOPEsw = SCOPEsw + scrolldata.VerticalScrollCount;
plotdatahelper();

function movePointers(varargin)
if get(gcf, 'userData') ~= -100
    pointerLoc = get(gcf, 'CurrentPoint');
    imageLoc = get(gca, 'Position');
    if length(get(gcf, 'userData')) > 1
        info = get(gcf, 'userData');
        blueRect = get(gcf, 'children');
        blueRect = blueRect(1);
        if info(1) == -1
            if pointerLoc(1) > info(2)
                set(blueRect, 'position', [info(2), .015, pointerLoc(1) - info(2), .01]);
            elseif pointerLoc(1) < info(2)
                set(blueRect, 'position', [pointerLoc(1), .015, info(2) - pointerLoc(1), .01]);
            end
        else
            if pointerLoc(2) > info(2)
                set(blueRect, 'position', [0, info(2), .01, pointerLoc(2) - info(2)]);
            elseif pointerLoc(2) < info(2)
                set(blueRect, 'position', [0, pointerLoc(2), .01, info(2) - pointerLoc(2)]);
            end
        end
    elseif pointerLoc(1) > imageLoc(1) && pointerLoc(1) < imageLoc(1) + imageLoc(3) && pointerLoc(2) > .05
        set(gcf, 'pointer', 'crosshair');
        kids = get(gcf, 'children');
        firstKid = get(kids(1), 'children');
        whichAxis = fix((pointerLoc(2) - .05) / .95 * (length(kids) - 2)) + 1;
        axisKids = get(kids(length(kids) - 1 - whichAxis), 'children');
        xCoord = (pointerLoc(1) - .05) / .95  * diff(get(gca,'Xlim')) + min(get(gca,'Xlim'));
        xData = get(axisKids(end - 4), 'xData');
        [junk whereX] = min(abs(xCoord - xData));

        % yCoord = min(abs(yCoords(xCoord) - round(((pointerLoc(2) - .05 - (whichAxis - 1) * (.95 / (length(kids) - 2))) / (.95 / (length(kids) - 2)))  * diff(get(kids(length(kids) - 1 - whichAxis),'ylim')) + min(get(kids(length(kids) - 1 - whichAxis),'ylim')))));

        if get(firstKid(length(firstKid)), 'userData')   == 0
            for index = 1:length(kids) - 2
                kidLines = get(kids(index), 'children');
                set(kidLines(length(kidLines)), 'xData', [xData(whereX) xData(whereX)], 'yData', get(kids(index), 'ylim'));
                yData = get(kidLines(length(kidLines) - 4), 'yData');
                yBounds = get(kids(index), 'ylim');
                if abs(yBounds(1) - yData(whereX)) > abs(yBounds(2) - yData(whereX))
                    whereY = yBounds(1) + .2 * diff(yBounds);
                else
                    whereY = yBounds(2) - .2 * diff(yBounds);
                end
                set(kidLines(length(kidLines) - 2), 'position', [xData(whereX) whereY], 'string', [' \bf' num2str(xData(whereX)) ', ' num2str(yData(whereX))]);
            end
        else
            for index = 1:length(kids) - 2
                kidLines = get(kids(index), 'children');
                set(kidLines(length(kidLines) - 1), 'xData', [xData(whereX) xData(whereX)], 'yData', get(kids(index), 'ylim'));
                yData = get(kidLines(length(kidLines) - 4), 'yData');
                firstPos = get(kidLines(length(kidLines) - 2), 'position');
                yBounds = get(kids(index), 'ylim');
                if abs(yBounds(1) - yData(whereX)) > abs(yBounds(2) - yData(whereX))
                    whereY = yBounds(1) + .1 * diff(yBounds);
                else
                    whereY = yBounds(2) - .1 * diff(yBounds);
                end
                set(kidLines(length(kidLines) - 3), 'position', [xData(whereX) whereY], 'string', [' \bf \Delta' num2str(xData(whereX)- firstPos(1)) ', \Delta' num2str(yData(whereX) - yData(1 + round(firstPos(1) / diff(xData(1:2)))))], 'visible', 'on');
            end
        end
    elseif pointerLoc(1) <= imageLoc(1) && pointerLoc(2) > 0.05
        set(gcf, 'pointer', 'top');
    elseif pointerLoc(1) > imageLoc(1) && pointerLoc(1) < imageLoc(1) + imageLoc(3) && pointerLoc(2) <= .05
        set(gcf, 'pointer', 'right');
    end
    if pointerLoc(2) < .015
        set(gcf, 'pointer', 'arrow');
    end
end

function mouseUpScope(varargin)
pointerLoc = get(gcf, 'CurrentPoint');
imageLoc = get(gca, 'Position');
if length(get(gcf, 'userData')) > 1
    info = get(gcf, 'userData');
    set(gcf, 'userData', 0);
    blueRect = get(gcf, 'children');
    delete(blueRect(1));
    kids = get(gcf, 'children');
    firstKid = get(kids(1), 'children');
    kidData = get(firstKid(end - 4), 'xData');
    if info(1) == -1
        if info(2) ~= pointerLoc(1)
            myPoint = (pointerLoc(1) - .05) / .95  * diff(get(kids(1),'Xlim')) + min(get(kids(1),'Xlim'));
            if myPoint < min(kidData)
                myPoint = min(kidData);
            end
            if myPoint > max(kidData)
                myPoint = max(kidData);
            end
            set(kids(1:length(kids) - 2), 'xlim', sort([(info(2) - .05) / .95  * diff(get(kids(1),'Xlim')) + min(get(kids(1),'Xlim')) myPoint]));
            set(kids(length(kids) - 1), 'string', num2str((max(kidData) - min(kidData)) / diff(get(kids(1), 'Xlim')),'%-1.1f'));
            xBounds = get(kids(1), 'xlim');

            zoomFactor = str2double(get(kids(length(kids) - 1), 'string'));
            newStep = 1 / zoomFactor / (1- 1 / zoomFactor);
            if newStep > 10
                set(kids(length(kids)), 'sliderStep', [1 newStep]);
            else
                set(kids(length(kids)), 'sliderStep', [newStep / 10 newStep]);
            end
            set(kids(length(kids)), 'value', xBounds(1) / max(kidData) / (1- 1 / zoomFactor));
        end
    else
        if info(2) ~= pointerLoc(2)
            set(kids(length(kids) - 1 - info(1)), 'ylim', sort((([info(2) pointerLoc(2)] - .05 - (info(1) - 1) * (.95 / (length(kids) - 2))) / (.95 / (length(kids) - 2)))  * diff(get(kids(length(kids) - 1 - info(1)),'ylim')) + min(get(kids(length(kids) - 1 - info(1)),'ylim'))));
        end
    end
elseif strcmp(get(gcf, 'SelectionType'), 'alt')
    kids = get(gcf, 'children');
    if pointerLoc(1) <= imageLoc(1) && pointerLoc(2) > 0.05
        whichAxis = fix((pointerLoc(2) - .05) / .95 * (length(kids) - 2)) + 1;
        lineKids = get(kids(length(kids) - 1 - whichAxis), 'children');
        tempXData = get(lineKids(1), 'xData');
        tempYData = get(lineKids(1), 'yData');
        tempBounds = get(kids(length(kids) - 1 - whichAxis), 'xlim');
        minBound = find(tempXData == min(tempBounds));
        maxBound = find(tempXData == max(tempBounds));
        set(lineKids(length(lineKids)), 'yData', [min(tempYData(minBound:maxBound)) max(tempYData(minBound:maxBound))]);
        set(lineKids(length(lineKids) - 1), 'yData', [min(tempYData(minBound:maxBound)) max(tempYData(minBound:maxBound))]);
        set(kids(length(kids) - 1 - whichAxis), 'YLimMode', 'auto');
    elseif pointerLoc(1) > imageLoc(1) && pointerLoc(1) < imageLoc(1) + imageLoc(3) && pointerLoc(2) <= .05
        set(kids(1:length(kids) - 2), 'XLimMode', 'auto');
        set(kids(length(kids) - 1), 'string', '1.0');
        set(kids(length(kids)), 'sliderStep', [1 1/0]);
    end
end

function traceScroll(varargin)

kids = get(gcf, 'children');
firstKid = get(kids(1), 'children');
xData = get(firstKid(end - 4), 'xdata');
zoomFactor = str2double(get(kids(length(kids) - 1), 'string'));
scrollValue = get(kids(length(kids)), 'value') * (1- 1 / zoomFactor);

windowSize = (xData(end) - xData(1)) / zoomFactor;
newStep = 1 / zoomFactor / (1- 1 / zoomFactor);
if newStep > 10
    set(kids(length(kids)), 'sliderStep', [1 newStep]);
else
    set(kids(length(kids)), 'sliderStep', [newStep / 10 newStep]);
end
set(kids(1:length(kids) - 2), 'Xlim', [xData(end) * scrollValue xData(end) * scrollValue + windowSize], 'dataaspectratiomode', 'auto', 'plotboxaspectratiomode', 'auto');

function    plotdatahelper()
global SCOPExData SCOPEdata SCOPEsw SCOPEhandles
persistent lhandle oldsw;

if isempty(oldsw);oldsw = -1;end;
if ndims(SCOPEdata)<3; SCOPEsw = 1; elseif SCOPEsw <1; SCOPEsw = 1; elseif SCOPEsw>size(SCOPEdata,3); SCOPEsw=size(SCOPEdata,3); end % check sweep limits
if oldsw ~= SCOPEsw % don't do anything if nothing has changed;
    if exist('lhandle'); if ~isempty(lhandle); delete(lhandle); end ;end
    for index = 1:size(SCOPEdata, 2); % nchns
        axes(SCOPEhandles(index));%cla;
        lhandle(index) = line(SCOPExData, SCOPEdata(:, size(SCOPEdata, 2) + 1 - index,SCOPEsw)); % plot SCOPEdata
        if max(SCOPEdata(:, size(SCOPEdata, 2) + 1 - index)) > min(SCOPEdata(:, size(SCOPEdata, 2) + 1 - index)) % set ylim
            set(gca, 'ylim', [min(SCOPEdata(:, size(SCOPEdata, 2) + 1 - index,SCOPEsw)) max(SCOPEdata(:, size(SCOPEdata, 2) + 1 - index,SCOPEsw))]);
        else
            set(gca, 'ylim', [SCOPEdata(1, size(SCOPEdata, 2) + 1 - index,SCOPEsw) - 1 SCOPEdata(1, size(data, 2) + 1 - index,SCOPEsw) + 1]);
        end
    end
end
oldsw = SCOPEsw;