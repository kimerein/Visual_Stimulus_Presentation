

clear R
for i=1:size(event,2)
    m= mean(event(i).isi(find(~isnan(event(i).isi))));
figure(99)
%% index of 
ind2 = find(event(i).isi>(m-10) &event(i).isi<(m+10) )
x = event(i).peak(ind2);
y = event(i).peak(ind2-1);
% plot(x,y,'.r')
conf = 0.01;
r= corrcoef(x,y,'alpha',conf);
R(i) = r(2,1);
title([' R:' num2str(R(i),'%1.2f') ' (p<' num2str(conf,'%1.2f') ')' ],'Interpreter','none');
end

i=3
    m= mean(event(i).isi(find(~isnan(event(i).isi))));
figure(99)
%% index of 
ind2 = find(event(i).isi>(m-10) &event(i).isi<(m+10) )
x = event(i).peak(ind2);
y = event(i).peak(ind2-1);
plot(x,y,'.r')
conf = 0.01;
r= corrcoef(x,y,'alpha',conf);
R(i) = r(2,1);
title([' R:' num2str(R(i),'%1.2f') ' (p<' num2str(conf,'%1.2f') ')' ],'Interpreter','none');
event(i).sfilename
