function h = fig(h,sname)
% function h = fig(h,sname)
% create new figure or call activate figure h and name it

if nargin <1 || isempty(h)
    h = figure;
else
    figure(h);
end
if nargin>2
    if ~isempty(sname)
        set(h,'name',sname);
    end
end