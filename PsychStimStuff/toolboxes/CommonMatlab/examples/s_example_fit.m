figure(100)
y = b(1:30);
x = [1:size(y)]';
p = polyfit(x,y,6);
xf = [0:0.1:size(x)]';
f = polyval(p,xf);
plot(x,y/1000,'ob',xf,f/1000,'-g')
hold on;
plot(xf(3:end),ddf,'-r')

figure(101)
plot(y)
hold all
plot(3:size(y,1),diff(diff(y)),'r')

%%% DOES NOT work don't know why
%% fits are variable
%% sum of Expontentials exp(x./-td)-exp(x./-tr);
% function fitdata = eventFit(y,riseWindow,fallWindow)
%
% const = cellstr(['tr';'td']);
ftype = fittype( 'exp(x/-td)-exp(x/-tr)','ind','x')
x = [riseWindow(1):riseWindow(2)]';
opts = fitoptions('Method','NonlinearLeastSquares','Algorithm','Gauss-Newton','Robust','on','StartPoint',[10 1])
[a gof out] = fit(x,-1*y(riseWindow(1):riseWindow(2)),ftype,opts)
% [a] = fit(x,-1*y(riseWindow(1):riseWindow(2)),ftype)

%% Naka rushton
% WHY do these two functions give very different fits?
Estimates=lsqcurvefit(@nakaRushtonFun,Starting,x,y)
[Estimates,R,J,COVB,MSE] =NLINFIT(x,y,@nakaRushtonFun,Starting)
Estimates = real(Estimates);
% To check the fit
fx = [min(x):min(x):max(x)];
figure(3);clf;
plot(x,y,'*')
hold on
plot(fx,Estimates(1).*(fx.^Estimates(3)./(Estimates(2).^Estimates(3) + fx.^Estimates(3))) + Estimates(4),'r')
