% mI = mNRespOdor
% sI = semNRespOdor
% mE = mNRespOdor
% sE = semNRespOdor
% 
save('cindy040608Neurodinner','mE','mI','sE','sI','con_epsc','con_ipsc','bac_epsc','bac_ipsc')
errorbar([1:size(mI,1)],mI(:,5),sI(:,5),'Color','b'); hold on;
errorbar([1:size(mE,1)],mE(:,5),sE(:,5),'Color','r')
ylim([0 1])
plotset(1)

gnames = {'0.25%','1%','5%','10%'};
figure(2);clf
barweb([mI(:,5)'; mE(:,5)']',[sI(:,5)' ;sE(:,5)']',1.75,gnames)
plotset(1)
% bar([mI(:,5)'; mE(:,5)']')
%%

figure(3);clf;
X = [1 2]
plot(ones(1,length(con_epsc))*X(1),con_epsc,'ro');hold on;
plot(ones(1,length(bac_epsc))*X(2),bac_epsc,'ro');hold on;
plot(X(1),mean(con_epsc),'k.','MarkerSize',30,'linewidth',3)
plot(X(2),mean(bac_epsc),'k.','MarkerSize',30,'linewidth',3)
plotJoinLine(con_epsc,bac_epsc,X,3);
xlim([0.5 2.5]);
set(gca,'XTick',[X])
set(gca,'XTickLabel',{'Ctrl';'Bac'})
plotset(1)

figure(4);clf;
X = [1 2]
plot(ones(1,length(con_ipsc))*X(1),con_ipsc,'bo');hold on;
plot(ones(1,length(bac_ipsc))*X(2),bac_ipsc,'bo');hold on;
plot(X(1),mean(con_ipsc),'k.','MarkerSize',30,'linewidth',3)
plot(X(2),mean(bac_ipsc),'k.','MarkerSize',30,'linewidth',3)
plotJoinLine(con_ipsc,bac_ipsc,X,3);
xlim([0.5 2.5]);
set(gca,'XTick',[X])
set(gca,'XTickLabel',{'Ctrl';'Bac'})
plotset(1)

