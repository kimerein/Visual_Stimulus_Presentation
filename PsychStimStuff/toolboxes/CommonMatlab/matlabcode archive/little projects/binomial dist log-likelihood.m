% find log-likelihood of binomial distribution
N = 14;
X = 7;
P = [0:0.01:1];
LL = -log(factorial(N)/(factorial(X)*factorial(N-X)).*P.^X.*(1-P).^(N-X));
plot(P,LL);
min_LL1 = min(LL)
maxP1 = P(find(LL==min_LL))
N = 9;
X = 8;
P = [0:0.01:1];
LL = -log(factorial(N)/(factorial(X)*factorial(N-X)).*P.^X.*(1-P).^(N-X));
hold on;
plot(P,LL,'r')
min_LL2 = min(LL)
maxP1 = P(find(LL==min_LL))
% http://fisher.forestry.uga.edu/popdyn/Likelihood.html

chi2inv(.95,2) % chisquared for 2 degree of freedom

