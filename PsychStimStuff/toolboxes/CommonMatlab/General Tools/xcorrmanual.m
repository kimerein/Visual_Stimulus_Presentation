%manual crosscorrelation
function [crosscorr indataX] = xcorrmanual(dataX,dataY,range,bnorm)
% function [crosscorr indataX] = xcorrmanual(dataX,dataY,range,bnorm)
%range should be an integer
if nargin<4
    bnorm = 0;
end
indataX{2*range +1}= []; %%
crosscorr = zeros(1,2*range +1);
for i = 1 :size(dataX,2)
    if dataX(i) % non zero
        for j = range:-1: -range
            jndex = i+j;
            if jndex > 0 && jndex <= size(dataY,2)
                temp = dataX(i)*dataY(jndex);
                crosscorr(j+range+1)  = crosscorr(j+range+1) + temp;
                if temp
                    indataX{j+range+1} = [indataX{j+range+1};i] ; %%  WARNING needs to be flipped in order like crosscorr on lin27
                end
            end
        end
    end
end

if bnorm %normalize so that correlation at zerolag is 1 (same as 'coeff' in xcorr
    crosscorr = crosscorr ./ (crosscorr(range+1));
end
crosscorr = crosscorr(end:-1:1);
% indataX{1:end} = indataX{end:-1:1};