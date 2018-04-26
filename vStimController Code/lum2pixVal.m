function p = lum2pixVal(l)
% function p = lum2pixVal(l)
% 
% INPUT
%   l: luminance
%   a: Scaling factor for power law relation
%   b: gamma exponent
%
% OUTPUT
%   p: pixel value
%

% Modified: SRO - 2/26/11

rigSpecific;

size(l)
a
b
p = round(exp(log(l./a)/b));     % a and b

if p > 255
    p = 225;
end

if p < 0 
    p = 0;
end