function revP =  nernstP(cIN,cOUT,z,tinC)
% function revP =  nernstP(cIN,cOUT,z,t)
if nargin <3 || isempty(z)
    z = 1; % charge
end
if nargin <4 || isempty(tinC)
    tinC = 30 ;
end

t = tinC  + 273.15; % Kelvin
R = 8.314472; % J · K-1 · mol-1
F = 96485; % J·V-1·mol-1
revP = - R*t/z/F*log(cIN/cOUT)*1000; %mV
