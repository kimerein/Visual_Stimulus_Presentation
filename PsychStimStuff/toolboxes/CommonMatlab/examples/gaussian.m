x = [-1000:1000];
s =50;
y = 1/(sqrt(2*pi)*s)*exp(-x.^2/s^2);
figure(1)
plot(x,y)

r = conv(y,y);
r2 = conv(r,r);
figure(2)
plot(r2)