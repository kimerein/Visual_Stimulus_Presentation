%% Sequence of neurons

%manual crosscorrelation
function [crosscorr] = seqmanual(dataX,dataY,dataZ,range)
%range should be an integer
if nargin<4
    bnorm = 0;
end

crosscorr = zeros(2*range +1,2*range +1);
for i = 1 :size(dataX,2)
    if dataX(i) % non zero
        for j = range:-1: -range
            jndex = i+j;
            if jndex > 0 && jndex <= size(dataY,2)
                for k = range:-1: -range
                    kndex = i+k;
                    if kndex > 0 && kndex <= size(dataZ,2)
                        temp = dataX(i)*dataY(jndex)*dataZ(kndex);
                        crosscorr(j+range+1,k+range+1)  = crosscorr(j+range+1,k+range+1) + temp;
                    end
                end
            end
        end
    end
end

crosscorr = crosscorr(end:-1:1);
crosscorr = reshape(crosscorr,2*range +1,2*range +1);