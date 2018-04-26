function f = Ferror(D)
global aaa;
global EP;
G = aaa'; %% GOAL ( compound EPSC)
%% NOTE conv(EP,D) and G must have same dimensions
%% EP globally defined (uEPSC)
c = conv(EP,D)';

% di = length(c) - length(G);
% if di >0
%     G =     
f = sqrt(sum((c - G).^2)); %% define error


