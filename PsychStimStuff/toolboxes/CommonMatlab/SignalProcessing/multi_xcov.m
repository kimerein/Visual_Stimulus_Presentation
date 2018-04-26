function covxy = multi_xcov(x,y,range)
% function covxy = multi_xcov(x,y,range)
% xcov on multiple sweeps (each row is a new sweep)
%   range:   to run xcov over (samples)
% BA 51508
range = round(range);
for i=1:size(x,1)
        covxy(i,:) =     xcov(x(i,:),y(i,:),range,'coeff')';
end