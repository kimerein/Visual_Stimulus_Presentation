function out =  alphan(x, Vrest)
% //an= -.0018 * (*x+50) / (exp (-(50+*x) / 25.) - 1);

out =0.032*(15-(x-Vrest))/(exp((15-(x-Vrest))/5)-1);

% //an=-0.01*(50+*x)/(exp(-(50+*x)/10)-1);
