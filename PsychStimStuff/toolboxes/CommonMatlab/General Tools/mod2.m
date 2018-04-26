function out = mod2(x,y)
% function out= mod2(x,y)
% BA function that is the same as mod
% but if x=y it returns y rather than 0

out = mod(x,y);
if out ==0
    out = y;
end