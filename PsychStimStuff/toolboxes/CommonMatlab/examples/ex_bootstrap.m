y = exprnd(5,100,1);
stats = bootstrp(100, @(x) [mean(x) std(x)], y);
plot(stats(:,1),stats(:,2),'o')