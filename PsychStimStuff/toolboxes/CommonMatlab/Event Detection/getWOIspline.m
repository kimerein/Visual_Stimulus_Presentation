function [subdata skipped dsx subdata_NOTspline]= getWOIspline(data,index,WOI,Nspline,bdetrend)
% function [subdata skipped dsx subdata_NOTspline]= getWOIspline(data,index,WOI,bdetrend))
%
% extract WOI from data around each of the index and spline to Nspline
% points
% INPUT: data  = data should contain 1 channel data from an experiment
%                each row should contain 1 sweep
%       in the same sweep)
%       index = array index in data (assume linear array)
%       WOI  = 1x2 array specifing window before after index to take
%       bdetrend = detrend each extracted window;
% OUTPUT: subdata = each row contains event around index
%         size(subdata) = [length(index) WOI(1)+WOI(2)+1 ]
%         skipped = vector of 1xlength(index)
%         if WOI does not exist for an index 1 will be placed in skipped.
%         dsx = sampling interval (relative to 1 sample in the input)
%                    i.e. dsx<1 if Nsamples > WOI(1)+WOI(2) +1
% BA103006
bdebug = 0;

if nargin <5;    bdetrend = 0; end

 data = data'; % switch so that each COL has one sweep (so indexes work)
if isrowvector(data); data = data'; end

x = [1:WOI(1)+WOI(2)+1];
dsx = x(end)/Nspline; sx =[dsx:dsx:x(end)]; 

skipped = zeros(1,length(index),'int16');
subdata =  zeros(length(index),Nspline);
if nargout > 3;    subdata_NOTspline =  zeros(length(index),sum(WOI)+1); end

ii=1;
% Ncycles = size(index,2);
Ncycles = length(index);
for i = 1: Ncycles
    try
        if floor((index(i) - WOI(1)-1)/size(data,1)) == floor((index(i) + WOI(2))/size(data,1)) % check in same sweep
            
            if bdetrend; temp = detrend(data(index(i)-WOI(1):index(i)+WOI(2)));                 
            else temp = data(index(i)-WOI(1):index(i)+WOI(2)); end

            if nargout > 3; subdata_NOTspline(ii,:) = temp; end
            subdata(ii,:)=pchip(x,temp,sx); % spline
            %        subdata(ii,:)=spline(x,temp,sx);

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

