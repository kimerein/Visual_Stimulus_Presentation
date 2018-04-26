function plotSMsem(data,xpos,pparam,bpsem)
% function plotSMsem(data,xpos,dt,pparam,bpsem)
% plot mean of data with sem, where data is a single value,
% INPUT: data - matrix, if it is a time series, time should run down
% columns
%        pparam - (opt) structure with plot parameters
%                  pparam.fid = <1>;
%                  pparam.color = <'b'>;
%        xpos - position of x axis
%        bpsem- <1> plot std if 1.
%BA031007

if nargin < 2
    xpos =1;
end
if nargin < 3 || isempty(pparam)
    fid = figure;
    pparam.color = 'k';
else
    if isfield(pparam,'fid')
        figure(pparam.fid);
    else
        fid = figure;
    end
end

if nargin < 4
    bpsem = 1;
end
linestyleM = ['s' pparam.color];

plot(xpos*ones(1,length(data)),data,['o' pparam.color],'MarkerSize',5 )
hold on
plot(xpos,mean(data),linestyleM)
hold on
if bpsem
    errorbar(xpos,mean(data),sem(data));
end
% axis tight