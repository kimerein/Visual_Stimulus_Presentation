function pparam = plotMstd(data,dt,pparam,bpstd,bMstd,bBase,bNorm)
% function pparam = plotMstd(data,dt,pparam,bpstd,bMstd,bBase,bNorm)
% plot mean of data with std
% INPUT: data - vector
%        dt  - scalar = 1/ sampling rate, or vector of times
%        pparam - (opt) structure with plot parameters
%                  pparam.fid = <1>;
%                  pparam.color = <'b'>;
%        bpstd - <1> plot std if 1.
%        bMstd - if 1 subtract mean std before plotting, <0> 
%        bNorm - if -1 normalize min of mean  to -1
%                if 1 normalize max of mean  to 1
%                default <0>
%        bBase - if 1 baseline to zero <0> 
%             
%        <0>
% BA031007
% Update BA051707 (added bMstd, bBaseNorm, non-scalar dt)
lw =2;
if nargin < 3 || isempty(pparam)
     pparam.fid = gcf;
     pparam.color = '';

else
    if isvar('pparam.linew'); lw = pparam.linew; end
    if isvar('pparam.fid'); figure(pparam.fid); else;      pparam.fid = gcf;  end
    if ~isvar('pparam.color');   pparam.color = '';  end
end
if nargin < 4 ;   bpstd = 1; end % plot std
if nargin < 5; bMstd = 0; end % if 1 subtract mean std before plotting
if nargin < 6 || isempty(bBase); bBase = 0; end % if 1 baseline to zero 
if nargin < 7 || isempty(bNorm); bNorm = 0; end 



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
scal =1; % normalization factor (defaults to 1)
if bBase
   data =     baselineCorrect(data);
end
if bNorm ==0
    if min(size(data))==1
        plot(xtime,data,linestyleM,'linewidth',lw)
    else
        plot(xtime,nanmean(data),linestyleM,'linewidth',lw)
    end
else
    if min(size(data))==1
        temp = data;
    else
        temp = nanmean(data);
    end
    [junk scal] = normalize(temp*bNorm);
    plot(xtime,temp*scal,linestyleM,'linewidth',lw)
% end
%     if min(size(data))==1
%         plot(xtime,bNorm*normalize(data*bNorm),linestyleM,'linewidth',lw)
%     else
%         temp = nanmean(data);
%         scal = 1./abs(max(temp));
%         plot(xtime,temp*scal,linestyleM,'linewidth',lw)
%     end
end

    hold on
if bpstd
    if bMstd
        % this code removes the average std before plotting
        plot(xtime,mean(data)+(std(data)-mean(std(data)))/scal,linestyleS)
        hold on
        plot(xtime,mean(data)-(std(data)-mean(std(data)))/scal,linestyleS)
    else
        plot(xtime,nanmean(data)+(nanstd(data))/scal,linestyleS)
        hold on
        plot(xtime,nanmean(data)-(nanstd(data))/scal,linestyleS)
    end
end
axis tight