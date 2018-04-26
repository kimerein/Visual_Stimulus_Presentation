%% Distribution of gamma LFP amplitudes
figure(20)
[a x] = hist(max2min,100);
stairs(x,a)
xedge = quadrent(a,x);
quad_ind = quandindex(max2min,xedge);

sel_ind = find(gammaClass(:,csubctr) & quad_ind(:,1));
tmSA = mean(eventtrigger(sel_ind,:),1);
tsSA = std(eventtrigger(sel_ind,:),1);
sel_ind = find(gammaClass(:,csubctr) & quad_ind(:,4));
tmLA = mean(eventtrigger(sel_ind,:),1);
tsLA = std(eventtrigger(sel_ind,:),1);
figure(21)
clf
plot(tmSA,'-b','linewidth',2)
hold on
plot(tmLA,'-r','linewidth',2)

figure(21); clf; 
figure(23); clf;
    figure(24); clf;
for i=1:4
    sel_ind = find(gammaClass(:,csubctr) & quad_ind(:,i));
    tmSA = mean(eventtrigger(sel_ind,:),1);
    tsSA = std(eventtrigger(sel_ind,:),1);
    rmSA(i) = mean(risetime(sel_ind));
    xtime = [1:length(tsSA)].*output.dt*1000;
    figure(21)
    plot(xtime,tmSA,'linewidth',2)
    hold all
    figure(23)
    plot(xtime,tsSA,'linewidth',2)
    hold all
    figure(24)
    plot(i,mean(Vmax(sel_ind)),'.b');
    hold on
    errorbar(i,mean(Vmax(sel_ind)),std(Vmax(sel_ind)));
end
axis tight
rmSA

