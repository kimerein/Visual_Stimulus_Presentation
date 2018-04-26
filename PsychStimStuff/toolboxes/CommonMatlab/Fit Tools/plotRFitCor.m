function [figHandle p r] = plotRFitCor(x,y,figNum)
%function [figHandle p r] = plotFitCor(x,y,figNum)
% robust fit and plot
% NOTE: confidence intervals don't work for this yet
% INPUT: x and y are vectors, figNum is scalar
%  
% BA 061407
N = length(x);
xf = [0:1/10*max(y)/abs(max(y)):max(abs(y))*max(x)/abs(max(y))]';

[b stats] = robustfit(x',y');
%% from tom Lane at mathworks
% Use robustfit to do least squares, get confidence intervals from nlparci
bint = nlparci(b,stats.resid,'cov' ,stats.robust_s,'alpha',.05);

% p = polyfit(x,y,1);
p = circshift(b,1);
f = polyval(p,xf);
figure(figNum)
clf
plot(xf,f,'-r')
% plot(b(1)+b(2)*xf,xf,'r-')
hold on;
plot(x,y,'.b');
conf = 0.05;
[r P] = corrcoef(x,y,'alpha',conf);
%   [R,P,RLO,RUP] = corrcoef(x,y,'alpha',conf);
title(['y = ' num2str(p(1),'%1.2f') 'x ' num2str(p(2),'%1.2f') '  ' num2str(bint(2,1),'%1.2f')  '<m<' num2str(bint(2,2),'%1.2f') '  R:' num2str(r(2,1),'%1.2f') ' (p<' num2str(P(2,1),'%1.3f') ') N=' num2str(N,'%d')],'Interpreter','none');
