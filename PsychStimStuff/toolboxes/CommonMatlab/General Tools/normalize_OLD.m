function [dataout scal] = normalize(datain)
%%EXTRA
%% Normalize Negative Peak to 1

if iscell(datain)
    for i =1 :size(datain,1)
        dataout{i,1} = datain{i,1};
        scal = -1./min(datain{i,1}(:,:,3));
        Scalar =  ones(size(datain{i,1},1),1)*scal;
        dataout{i,1}(:,:,3) = datain{i,1}(:,:,3).*Scalar ;
    end
else
    dataout = datain;
    scal = -1./min(datain(:,:,3));
    Scalar =  ones(size(datain,1),1)*scal;
    dataout(:,:,3) = datain(:,:,3).*Scalar ;
end