function av = getfig(handles)
% function av = getfig(handles)
%
% gets handles for figure object (copied from figs2subplots.m)
%
% BA071507

% Parsing handles vector
av = [];
for k = 1:length(handles)
    if strcmp(get(handles(k),'Type'),'axes')
        av = [av handles(k)];
    elseif strcmp(get(handles(k),'Type'),'figure');
        fc = get(handles(k),'Children');
        for j = length(fc):-1:1
            if strcmp(get(fc(j),'Type'),'axes') && ~strcmp(get(fc(j),'Tag'),'legend')
                av = [av fc(j)];
            end;
        end;
    end;
end;