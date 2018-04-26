function mass = getmass(molarmass,con,vol)
% function mass = getmass(molarmass,con,vol)
% INPUT:
%          molarmass (g/mol)
%          con (mol)
%          vol (liters)
%
% mass from concentration and volume
% Solves for mass
% mass (g)
% ---                  * 1/Volume = Concentration
% molar_mass (g/m)

mass = molarmass*con*vol;

