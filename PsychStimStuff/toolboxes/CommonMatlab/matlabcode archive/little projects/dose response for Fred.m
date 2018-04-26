% for FP
% dose response curve
A1 = 0;
A2 = 100;
x0 = 49.981;
x = [-100:1:100];
p = .04403;
D = A1 + (A2 -A1)./(1 + 10.^(((x0) - x).*p));

figure(99);
plot(x,D)
hold all;