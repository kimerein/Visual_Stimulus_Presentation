function out = betam(x,Vrest)
% //bm=.12 * (*x - 5.9) / (exp ((*x - 5.9) / 4.47) - 1);

out=.28*((x-Vrest)-40)/(exp(((x-Vrest)-40)/5)-1);
% /*bm=4*exp(-(60+*x)/18);*/
