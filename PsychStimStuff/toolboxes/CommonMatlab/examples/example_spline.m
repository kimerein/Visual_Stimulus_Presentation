figure(100)
y = mean(EE(5:100,:),2);
x = [1:size(y)]';
step = 1; xs = [1:step:size(y)]';
s = spline(x,y,xs);
p = pchip(x,y,xs);
plot(x,y,'-b',xs,s,'-g',xs,p,'-r')
plot(x(1:end-1),y(1:end-1),'-b',x(1:end-1),diff(y),'-g')

% hold on;
% plot(xf(3:end),ddf,'-r')
% 
% figure(101)
% plot(y)
% hold all
% plot(3:size(y,1),diff(diff(y)),'r')