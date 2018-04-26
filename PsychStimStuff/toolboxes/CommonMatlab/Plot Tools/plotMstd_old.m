function pparam = plotMstd(data,dt,pparam,bpstd)
% function plotMstd(data,dt,pparam,bpstd)
% plot mean of data with std
% INPUT: data - vector
%        dt  - 1/ sampling rate
%        pparam - (opt) structure with plot parameters
%                  pparam.fid = <1>;
%                  pparam.color = <'b'>;
%        bpstd - <1> plot std if 1.
%BA031007
lw =2;
if nargin < 3 || isempty(pparam)
    pparam.fid = figure();
    pparam.color = 'b';
        
else
    figure(pparam.fid);
        if isvar('pparam.linew'); lw = pparam.linew; end
end
if nargin < 4
    bpstd = 1;
end
linestyleM = ['-' pparam.color];
linestyleS = ['--' pparam.color];

if max(size(dt))==1
xtime = [1:size(data,2)]*dt*1000; % converts to ms from s
else
    if length(dt) ~= size(data,2)
        error('dt must a scalar, or a vector whose length is equal to the number of columns in "data"');
    end
    xtime = dt;
end
% toc
if min(size(data))==1
    plot(xtime,data,linestyleM,'linewidth',lw)
else
    plot(xtime,nanmean(data),linestyleM,'linewidth',lw)
end
hold on
if bpstd
    % this code removes the average std before plotting
    plot(xtime,mean(data)+(std(data)-mean(std(data))),linestyleS)
    hold on
    plot(xtime,mean(data)-(std(data)-mean(std(data))),linestyleS)
%     plot(xtime,nanmean(data)+(nanstd(data)),linestyleS)
%     hold on
%     plot(xtime,nanmean(data)-(nanstd(data)),linestyleS)
end
axis tight