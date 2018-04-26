function stats = expFit(y,dt,win)
% function stats = expFit(y,dt)
% dt - (optional)
% win = (optional) window for fit (in samples)
bdebug =1;
if nargin<2 || isempty(dt)
    dt = 1;
end
if nargin<3 || isempty(win)
    win = [1 length(y)];
end
if isrowvector(y)
    y = y';
end
x = [win(1):win(2)]*dt;
y = y-min(min(y))+1; % insure all y is positive
ly = log(y);
for i=1:size(ly,2)
%     [b,bint,r,rint,ostats] 
    [b bint] = regress(ly(win(1):win(2),i),[ones(length(x),1) x']);
    % p = polyfit(x',y(:,i),1);
    stats(i).p = circshift(b,1)';
    stats(i).pconf = bint;
end
if bdebug
    xtime  = [1:length(y)]*dt;

    for i=1:size(ly,2)
%         figure(3000);clf;
%         semilogy(xtime,y(:,i)); hold on;
%         semilogy(x,exp( stats(i).p(2)).*exp( stats(i).p(1)*x),'-k','linewidth',2)
%         axis tight

        figure(3001);
        plot(xtime,y(:,i)); hold on;
        plot(x,exp( stats(i).p(2)).*exp( stats(i).p(1)*x),'-r','linewidth',2)
        axis tight
%         pause
    end
end

