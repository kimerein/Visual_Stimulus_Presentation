function [alpha k] = alphafit(x,y)
% function [alpha k] = alphafit(x,y)
% fits x,y row/column vectors with alpha function
% and returns 

%% make col vects into row vects for fit()
if size(y,2) ~= 1
    y = y';
end
if size(x,2) ~= 1
    x = x';
end

ftype = fittype( 'k*x*exp(-alpha*x)','ind','x')
% opts = fitoptions('Method','NonlinearLeastSquares','Algorithm','Gauss-Newton','Robust','on','StartPoint',[10 1])
% [a gof out] = fit(x,y,ftype,opts)
[a gof out] = fit(x,y,ftype);
temp = coeffvalues(a);
alpha = temp(1);
k = temp(2);