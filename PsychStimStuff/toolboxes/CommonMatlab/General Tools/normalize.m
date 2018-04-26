function [dataout scal] = normalize(datain, range)
% function [dataout scal] = normalize(datain, range)
% normalizes  datain to negative peak 
%    datain should have a single sweep on each row
% range is optional arguement containing 2 element vector with range to be
% used for normalized
%  e.g. range =  [5 300], data would be baselined to the mean value of each row between indices of
%  5 and 300
% BA051508

if nargin <2  % default to baseline whole sweep
    range = [1 size(datain,2)];
end

if size(datain,2)==1 % turn into rowvector
    datain = datain';
end
if length(range)==1
    scal = 1./abs(datain(:,range));
else
    scal = 1./abs(min(datain(:,range(1):range(2))')');
end
dataout =  datain.*repmat(scal,1,size(datain,2));

%% OLD VERSION FOR something (don't remember what
%% Normalize Negative Peak to 1
% 
% if iscell(datain)
%     for i =1 :size(datain,1)
%         dataout{i,1} = datain{i,1};
%         scal = -1./min(datain{i,1}(:,:,3));
%         Scalar =  ones(size(datain{i,1},1),1)*scal;
%         dataout{i,1}(:,:,3) = datain{i,1}(:,:,3).*Scalar ;
%     end
% else
%     dataout = datain;
%     scal = -1./min(datain(:,:,3));
%     Scalar =  ones(size(datain,1),1)*scal;
%     dataout(:,:,3) = datain(:,:,3).*Scalar ;
% end