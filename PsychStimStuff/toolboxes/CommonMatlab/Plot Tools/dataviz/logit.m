function y = logit(x)
%  apply logit transform
%  y = logit(x)
%  where logit is defined by y = log(x/(1-x))
%  intended for 0<x<1

% Copyright (c) 1998-2000 by Datatool
% $Revision: 1.10 $

%  deal with negative x if necessary
x = max(x,0);

%  deal with large x if necessary
y = (1/eps)*ones(size(x));
y(x<1) = x(x<1)./(1-x(x<1));

%  deal with nonpositive y if necessary
y(y<=0) = eps;

y = log(y);
