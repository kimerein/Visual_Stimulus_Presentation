function vol = getvol(molarmass,con,mass)
% function mass = getmass(molarmass,con,vol)
% INPUT:
%          molarmass (g/mol)
%          con (mol)
%          mass (g)
%
% vol from concentration and mass
% Solves for mass
% vol(L)
% ---                  * 1/Volume = Concentration
% molar_mass (g/m)
 vol = mass/molarmass/con;
% mass = molarmass*con*vol;
