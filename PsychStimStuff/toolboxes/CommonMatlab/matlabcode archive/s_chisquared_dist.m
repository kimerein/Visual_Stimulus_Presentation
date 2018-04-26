%% 
demonstration of chi-squared distribution
[a x] = hist(randn(1e5,1).^2 + randn(1e5,1).^2+ randn(1e5,1).^2+ randn(1e5,1).^2+ randn(1e5,1).^2,100);
figure(99);clf;
plot(x,a/max(a))
hold on
x2 = [0:.01:15];
y2 = chi2pdf(x2,5);
plot(x2,y2/max(y2))