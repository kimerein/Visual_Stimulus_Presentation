function [ind dist] = match(x,y,w)
%function ind = match(x,y)
% x is a scalar, y is a vector
% function finds the value index of an element in y that is within w of x
%% if no such value is found -1 is returned
temp = min(abs(y - x));
temp2 = abs(y-x);

if temp <= w
    ind = find(temp2 == temp,1,'first');
else 
    ind = -1;
end

dist = temp;
