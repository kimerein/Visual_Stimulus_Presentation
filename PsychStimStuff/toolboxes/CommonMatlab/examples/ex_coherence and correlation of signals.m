cyc =10;
x = [0:0.01:cyc*2*pi];
% rn times a cycle change phase
rn = 2;
rphase = rand(1,cyc*rn)*2*pi;
r = [];
for j = 1:length(rphase)
    r = [r rphase(j)*ones(1,length(x)/cyc/rn)];
end

x = x(1:length(r));
s = sin(x);
s1 = sin(x+pi/3)*sin(5.3*x);
s2 = sin(x +r);
figure(298);clf;
plot(s);hold on ;plot(s1,'-k')%plot(s2,'-r');
figure(299);clf;
plot(xcorr(s2),'-k'); hold on;
plot(xcorr(s),'-k')
plot(xcorr(s2,s),'-g')

figure(301);clf;
plotcoh(s,s2,1,301,[],0,[]);
xlim([0 0.5])

figure(302);clf;
plotcoh(s,s1,1,302,[],0);
xlim([0 0.5])

s3 = sin(3*x);
figure(301);clf;
plot(xcorr(s3),'-k'); hold on;
plot(xcorr(s),'-k')
plot(xcorr(s3,s),'-g')

plotxcov(s,s,length(x)/cyc,1,299);
figure(300);
plotcoh(x,r,1,300);

figure;

r = rand(1,length(x))+10;
r = ones(1,length(x));

c = xcov(x,r);
figure(99)
plot(c)



%% TEST WITH coherent example
% WHY such poor coherence?
dX = 0.05;
N = 100;
X = [0:dX:N*2*pi];
S = sin(X)+ sin(1.4*X)+ sin(13*X); +5*randn(size(X));
P = [1:N]*2*pi;
% P =  P + randn(size(P))*2*pi/8;
P = round(P/dX);
event.spikes = P;
[C,phi,S12,S1,S2,f] = coherencycpt(S',P') ;

figure(1)
subplot(2,1,1)
plot(f,C);
subplot(2,1,2)
plot(f,phi);

