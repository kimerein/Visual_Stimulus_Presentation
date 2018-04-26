function [dataout baseadj] = baselineCorrect(datain,range)
% function [dataout baseadj] = baselineCorrect(datain,range)
% datain should have a single sweep on each row
% range is optional arguement containing 2 element vector with range to be
% used for baselining
%  e.g. range =  [5 300], data would be baselined to the mean value of each row between indices of
%  5 and 300
% BA051508

if size(datain,2)==1 % turn into rowvector
    datain = datain';
end
if nargin <2  % default to baseline whole sweep
    range = [1 size(datain,2)];
end
range = round(range);
if length(range)==1
    range = [range range];
end
%% shift baseline to zero
baseadj = mean(datain(:,range(1):range(2)),2);
dataout =  datain - repmat(baseadj,1,size(datain,2));


%% OLD VERSION FOR something (don't remember what
% if iscell(datain)
%     for i =1 :size(datain,1)
%         dataout{i,1} = datain{i,1};
%         baseadj = mean(datain{i,1}(1:int16(size(datain{i,1},1)*.1),:,3));
%         dataout{i,1}(:,:,3) =  datain{i,1}(:,:,3) - ones(size(datain{i,1}(:,:,3),1),1)*baseadj;
%     end
% else
% dataout = datain;
% baseadj = mean(datain(1:int16(size(datain,1)*.1),:,3));
% dataout(:,:,3) =  datain(:,:,3) - ones(size(datain(:,:,3),1),1)*baseadj;
% end