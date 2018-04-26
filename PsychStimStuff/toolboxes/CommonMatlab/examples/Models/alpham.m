function out = alpham(x, Vrest)
% //am=-1.74 * (*x - 11.0) / (exp (-(*x - 11.) / 12.94) - 1);

out=.32*(13-(x-Vrest))/(exp((13-(x-Vrest))/4.)-1); %  /*Top from RA neuron model used in archi paper */
