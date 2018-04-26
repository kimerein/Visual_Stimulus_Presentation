function out =  betan(x, Vrest)
% //bn=.0036 * (*x - 10) / (exp ((*x - 10) / 12) - 1);

out=.5*exp((10-(x-Vrest))/40);
% /*bn=.125*exp(-(60+*x)/80);*/
