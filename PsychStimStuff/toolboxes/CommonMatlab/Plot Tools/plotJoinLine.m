function fid = plotJoinLine(data1,data2,X,fid,brow)
% function fid = plotJoinLine(data1,data2,X,fid,brow)
%
% if brow = 1 then data will be plotted in row (default = 0)
%
% BA070708
if length(data1) ~= length(data2)
    error('length(data1) ~= length(data2');
end
if nargin<3
    X = [1 2];
end
if nargin<4
    fid = figure;
end
if nargin < 5
    brow = 0;
end

for i=1:length(data1)
    hold on;
    if ~brow
        line([X(1),X(2)],[data1(i),data2(i)],'Color','k');
    else
        line([data1(i),data2(i)],[X(1),X(2)],'Color','k');
    end
end
