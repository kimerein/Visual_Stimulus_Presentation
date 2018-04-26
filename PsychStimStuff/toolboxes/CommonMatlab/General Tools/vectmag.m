function out = vectmag(V)
% function out = vectmag(V)
% find length of vector
% out = sqrt(sum(V.^2))
S= size(V);
if length(S) >2;
    error(' vectmag function only takes 1-D vectors')
end
out = sqrt(sum(V.^2));
    