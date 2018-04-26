function out = gammaf(VV,theta,sigma)
% function out = gammaf(VV,theta,sigma)
out =1.0/(1.0+exp(-(VV-theta)/sigma));