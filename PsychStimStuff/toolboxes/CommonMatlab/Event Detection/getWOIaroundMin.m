function [subdata skipped index_min]= getWOIaroundMin(data,index,WOI)
% function [subdata skipped index_min]= getWOIaroundMin(data,index,WOI)
%
% extract WOI from data around each of the minimum. The minimum is detected within a window WOI of the
% index supplied.
%
% USAGE: to extract spikes centered at their minimum

% INPUT: data  = data should contain 1 channel data from an experiment
%                each row should contain 1 sweep
%       in the same sweep)
%       index = array index in data (assume linear array)
%       WOI  = 1x2 array specifing window before after index to take
% OUTPUT: subdata = each row contains event around index
%         size(subdata) = [length(index) WOI(1)+WOI(2)+1 ]
%         skipped = vector of 1xlength(index)
%         if WOI does not exist for an index 1 will be placed in skipped.
%         index_min = of minimum
%                         note: index_min may NOT be the same size as index
%                         if some indices are skipped
% BA012109 based on getWOI
%
bdebug = 0;
data = data'; % switch so that each COL has one sweep (so indexes work)
if isrowvector(data); data = data'; end

skipped = zeros(1,length(index),'int16');
subdata =  zeros(length(index),sum(WOI)+1);
index_min = zeros(size(index));
ii=1;
% Ncycles = size(index,2);
Ncycles = length(index);
for i = 1: Ncycles
    try
        if floor((index(i) - WOI(1)-1)/size(data,1)) == floor((index(i) + WOI(2))/size(data,1)) % check in same sweep

            temp = data(index(i)-WOI(1):index(i)+WOI(2));
            temp2 = find(temp==min(temp),1,'last');
            index_min(i) = temp2(end) + index(i) -(WOI(1)-1);

            subdata(ii,:) = data(index_min(i)-WOI(1):index_min(i)+WOI(2));
            % NOTE: not catching case where WOI doesn't exist around index_mind
            ii=ii+1;

            if bdebug
                figure(7)
                clf
                plot(data(index(i)-WOI(1):index(i)+WOI(2)))
            end
        else
            skipped(i) = 1;
        end
        if mod(i,50) ==0 %% Print progress
            sprintf([num2str(i) ' of ' num2str(Ncycles)]);
        end
    catch
        sprintf(['ERROR at i=' num2str(i,'%d') ' of ' num2str(Ncycles)])
        keyboard
    end

end

%% get rid of skippe
s = sum(skipped);
subdata = subdata(1:end-s,:);
index_min = index_min(1:end-s);

