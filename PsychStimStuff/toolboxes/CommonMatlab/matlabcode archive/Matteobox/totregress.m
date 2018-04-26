function b = totregress(y,x,dimout)
% totregress total least squares linear regression on 1-dimensional data
% (i.e., linear regression based on orthogonal distances)
%
% y = b(1) + b(2)*x
%
% b = totregress(y,x,dimout)

% Make sure all inputs are given
if nargin < 3
   dimout = [1 2];
end

% Make sure the data comes in row vectors
if size(x,1) < size(x,2)
   x = x';
end
if size(y,1) < size(y,2)
   y = y';
end

% The data
M = [x y];

% The means
m = mean(M)';

% The zero-mean data
N = M - repmat(m',[length(x) 1]);

% Compute the singular value decomposition
[U,S,V] = svd(N);

% The direction orthogonal to the best linear fit
v = V(:,2);

% The offset along the orthogonal direction
u0 = dot(m,v);

% The direction of maximal variance
p = [0 -1; 1 0]*v;

% Get the regression line
if p(1) == 0;
   disp('Degenerate case, is vertical line');
   B = [v(1)*u0 NaN];
else
   B_slope = p(2)/p(1);
   B_inter = v(2)*u0 - v(1)*u0*p(2)/p(1);
   B = [B_inter B_slope];
end

% Keep the outputs that are needed
b = B(dimout);

return

%% A simple example

% Intercept and slope of a line
k = [100 0.5];

% A noiseless line
x0 = 1:100;
y0 = k(1) + k(2)*x0;

% Add noise
x = x0 + 30*(rand(size(x0))-0.5);
y = y0 + 30*(rand(size(y0))-0.5);

% Total linear regression
b = totregress(y,x);

% Plot it
xlim = [min(x),max(x)];
figure; hold on;
plot(x,y,'bo');
plot(xlim,b(1)+b(2)*xlim,'r-');


%--- The relation between usual and total linear regression ---
ns = 1000;
vm = 100;
v1 = linspace(-vm,vm,ns) + (rand(1,ns)-0.5)*vm;
v2 = linspace(-vm,vm,ns) + (rand(1,ns)-0.5)*vm;
bv = regress(v2',[ones(ns,1), v1'], 0.05);
cv = totregress(v2,v1);
m = 0.5;
k = 2;
r1 = vm*(1-exp(-(max(0,v1)/(vm*m)).^k));
r2 = vm*(1-exp(-(max(0,v2)/(vm*m)).^k));
br = regress(r2',[ones(ns,1), r1'], 0.05);
cr = totregress(r2,r1);

figure; ax = [];
ax(1) = subplot(1,2,1); hold on; p = [];
axlimv = [min([v1 v2]),max([v1 v2])];
plot(v1,v2,'ko');
p(1) = plot(axlimv,bv(1)+bv(2)*axlimv,'r-');
p(2) = plot(axlimv,cv(1)+cv(2)*axlimv,'b-');
set(p,'linewidth',2);
set(gca,'xlim',axlimv,'ylim',axlimv,'plotbox',[1 1 1]);

ax(2) = subplot(1,2,2); hold on; p = [];
axlimr = [min([r1 r2]),max([r1 r2])];
plot(r1,r2,'ko');
p(1) = plot(axlimr,br(1)+br(2)*axlimr,'r-');
p(2) = plot(axlimr,cr(1)+cr(2)*axlimr,'b-');
set(p,'linewidth',2);
set(gca,'xlim',axlimr,'ylim',axlimr,'plotbox',[1 1 1]);

