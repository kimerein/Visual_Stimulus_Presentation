% example
% Image histogram

figure;
temp = imdf(:,:,1);
hist(double(temp(:)),1000)