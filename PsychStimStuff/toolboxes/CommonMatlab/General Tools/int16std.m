function std = int16std(data)
%function std = int16std(data)
% 
%   data can have up to 2 dimensions
%   
% Computes the standard deviation of integer data
% on dim 1
%

tsum = 0;
 temp1 = single(mean(data,1)); %mean 
for( i = 1:size(data,2)) %note this can be done with out for loops but require a large array of singles for the sqrt
    x = temp1(i);
%     for( ii = 1:size(data,1))
%         temp = (single(data(ii,i)-x))^2;
%         tsum = [tsum + temp];
%         if int32(ii/10000) ==1
%             ii
%         end
%     end
     temp = (single(data(:,i))-x.*ones(size(data,1),1,'single')).^2;
     tsum= double(sum(temp));
    std(i) = sqrt(tsum./((size(data,1)-1)));
end
