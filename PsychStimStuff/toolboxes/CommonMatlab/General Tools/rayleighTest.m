function [P Z]= rayleighTest(rdata,n)
% rayleghTest(rdata,n)
% rdata vector 
% n =  number different points on circle used to compute rdata %** do they need to be spaning an entire circle?
% if null hypothesis is rejected, we may conclude that there is a mean population direction, and if not, we may conclude that population distribution to be uniform around the circle;    
%  Assumption:  The population does not have more than one mode.
%
% BA 110309

% compute Rayleigh Z
Z = n*rdata.^2;

% P of R ( don't know if this is an approxmation from p70 "Statistical
% Analysis of Circular DAta N.I Fisher

P = exp(-Z).*(1+(2*Z-Z.^2)./(4*n)-(24*Z-132*Z.^2  + 76*Z.^3-9*Z.^4)/(288*n^2));