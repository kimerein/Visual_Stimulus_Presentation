function [ind_xthres sw] = thresh(data, threshold,minThres)
% function [ind_xthres] = thresh(data, threshold,minThres)
% finds point of threshold crossing (negative going)
%ind_xthres = find(diff(data < thres)==1)+1;
% minThres (opt) - min number of samples that must be above threshold to
% count as threshold crossing
% BA

if nargin >2 & ~isempty(minThres) % impose min number of samples that must be above threshold for threshold crossing to be included
    a = int16(data < threshold);
     for i=1:minThres
         a = a+ [a(i+1:end,:);  zeros(i,size(a,2),'int16')];
     end
     data(a<=minThres) = mean(data(:)); % remove events that don't meet minThres criterea
     clear a;
end

if nargout <2
    ind_xthres = find(diff(data < threshold)==1)+1; 
else
    [ind_xthres sw] = find(diff(data < threshold)==1);
    ind_xthres = ind_xthres +1;   
end

