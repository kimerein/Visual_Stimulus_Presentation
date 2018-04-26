function [amp instf] = philbert(data,dt)
% function [amp instf] = philbert(data,dt)
%plots phase angle and magnitude of analytic signal
if nargin <2
    dt = 1;
end
% A = hilbert(data,20000);
A = hilbert(data);

x = [1:length(data)]*dt;
figure(1)
subplot(2,2,1)
plot(x,data)
hold all
title('signal')
axis tight

subplot(2,2,2)
plot(A)
hold all
title('analytic')
axis tight

subplot(2,2,3)
amp = abs(A);
plot(x,amp)
hold all
title('magnitude')
axis tight

subplot(2,2,4)
instf = diff(unwrap(angle(A)))./dt;%./(2*pi);
plot(x(1:end-1),instf)
title('inst freq')
hold all
axis tight


