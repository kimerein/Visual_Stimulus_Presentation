function [Sout t]= reduceDim(S,t)
% function [S t]= reduceDim(S,t)
% converts 3 dimensional matrix S (e.g. time x frequency x trial) 
% into time x freq where each trials are concatinated in time.
% OR  2 dim  matrix to (e.g. time x sweeps) converts to vector
%
% USAGE: used in combitionat with chronux package mtspecgramc.m to convert
% data into form compatible with plotmatrix.m
% BA053008
if nargin <2; t = [1 2]; end
   
if length(size(S))==3
    Sout = zeros(size(S,1)*size(S,3),size(S,2));
    for i =1:size(S,3)
        Sout((i-1)*size(S,1)+1:(i-1)*size(S,1)+size(S,1),:) = S(:,:,i);
    end
elseif     length(size(S))==2
      Sout = S(1:end);
end

    t = [0:size(Sout,1)-1]*(t(2)-t(1))+t(1); % create time vector of correct length
