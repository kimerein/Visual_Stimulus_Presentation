ii=1; jj=1; clear exr; clear ixr;
for i=1:size(event,2)
    temp = find(~isnan(event(i).isi))
    x = event(i).isi(temp);
    y=circshift(event(i).peak(temp),1);
    conf = 0.01;
    N=length(x);
    p = polyfit(x,y,1);
    xf = [0:1/10*max(x)/abs(max(x)):max(abs(x))*max(x)/abs(max(x))]';
    f = polyval(p,xf);;
    figure(1)
    clf
    plot(xf,f,'-b')
    hold on;
    plot(x,y,'.b');
title(['y = ' num2str(p(1),'%1.2f') 'x ' num2str(p(2),'%1.2f') ' R:' num2str(r(2,1),'%1.2f') ' (p<' num2str(conf,'%1.2f') ') N=' num2str(N,'%d')],'Interpreter','none');
    if event(i).V < -50
        r = corrcoef(x,y,'alpha',conf);
        exr(ii) = r(2,1);
        ii = ii +1;
    else
        r = corrcoef(x,y,'alpha',conf);
        ixr(jj) = r(2,1);
        jj=jj+1;
    end
    pause;

end
exr
ixr
j= j+1
mxr(j,:) = [mean(exr) mean(ixr)] 