function [subdata skipped]= getWOI(data,index,WOI)
% function [subdata skipped]= getWOI(data,index,WOI)
%
% extract WOI from data around each of the index
% INPUT: data  = data should contain 1 channel data from an experiment
%                each row should contain 1 sweep
%       in the same sweep)
%       index = array index in data (assume linear array)
%       WOI  = 1x2 array specifing window before after index to take
% OUTPUT: subdata = each row contains event around index
%         size(subdata) = [length(index) WOI(1)+WOI(2)+1 ]
%         skipped = vector of 1xlength(index)
%         if WOI does not exist for an index 1 will be placed in skipped.
%
% BA103006
bdebug = 0;
 data = data'; % switch so that each COL has one sweep (so indexes work)
if isrowvector(data); data = data'; end

skipped = zeros(1,length(index),'int16');
subdata =  zeros(length(index),sum(WOI)+1);
ii=1;
% Ncycles = size(index,2);
Ncycles = length(index);
for i = 1: Ncycles
    try
        if floor((index(i) - WOI(1)-1)/size(data,1)) == floor((index(i) + WOI(2))/size(data,1)) % check in same sweep

            subdata(ii,:) = single(data(index(i)-WOI(1):index(i)+WOI(2)));
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

