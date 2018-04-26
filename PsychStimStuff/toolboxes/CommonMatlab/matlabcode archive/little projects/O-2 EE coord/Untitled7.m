close 10
figure(10)
plot(proData{4}(:,3))
hold on;
plot(smooth(proData{4}(:,3),18,'rloess'),'-r')
% plot(smooth(proData{5}(:,3),7,'sgolay'),'-r')
figure(101)
tic
y = smooth(proData{3}(:,3),15,'rloess');
x = [1:size(y)]';
step = .1; xs = [1:step:size(y)]';
s = spline(x,y,xs);
p = pchip(x,y,xs);
plot(x,proData{3}(:,3),'-b',xs,s,'-g');
toc

%% Extract spikes
%% Extract Parameters
%% data read in with "plotSelectedData.m"
for i = 1:size(proData,1)
    a = proData{i,1}(:,1,3);
    FSparam(i,:) = xSpikeParams(a');
end

%% Extract spikes
%% Extract Parameters
%% data read in with "plotSelectedData.m"
for i = 1:size(proData,1)
    a = proData{i,1}(:,1,3);
    PYparam(i,:) = xSpikeParams(a');
 end
%         [ 1  2  3  4   5   6   7   8   9  10]
%% parm = [P1 T1 P2 P1W T1W P2W T1F T1R P2R P2F];
close 100;
figure(107);
hold on;
for i = 1:size(FSparam,1)
%           plot((FSparam(i,8)+FSparam(i,6))*.002,abs(FSparam(i,2)/FSparam(i,3)),'.b')  %102
          plot((FSparam(i,6))*.002,abs(FSparam(i,2)/FSparam(i,3)),'.b')  %102
end
for i = 1:size(PYparam,1)
%           plot((PYparam(i,8)+PYparam(i,6))*.002,abs(PYparam(i,2)/PYparam(i,3)),'.r')  %102
          plot((PYparam(i,6))*.002,abs(PYparam(i,2)/PYparam(i,3)),'.r')  %102
end
figure(106)
%           plot((FSparam(:,8)+FSparam(:,6)).*.002.*abs(FSparam(:,2)./FSparam(:,3)),1,'.b',(PYparam(:,8)+PYparam(:,6)).*.002.*abs(PYparam(:,2)./PYparam(:,3)),1,'.r')  %102
%           plot((FSparam(:,8)+FSparam(:,6)).*.002,0,'.b',(PYparam(:,8)+PYparam(:,6)).*.002,0,'.r','MarkerSize',20)  %102
          plot((FSparam(:,5)).*.002,0,'.b',(PYparam(:,10)).*.002,0,'.r','MarkerSize',20)  %102
          ylim([-0.2 0.2])

params = [FSparam;PYparam];
figure(107)
          plot((params(:,8)+params(:,6)).*.002,abs(params(:,2)./params(:,3)),'.b','MarkerSize',20)  %102
figure(108)
          plot((params(:,8)+params(:,6)).*.002,ones(size(params,1),1),'.b','MarkerSize',20)  %102
figure(109)
          plot((params(:,8)).*.002,zeros(size(params,1),1),'.b','MarkerSize',20)  %102

figure(108)
          plot(b,)  %102

% if pca_flag ==1
[PC, SCORE, LATENT, TSQUARE] = princomp(single(params));
% end
% Scatter plot first two PCs
figure(10)
hold off;
set(gcf,'Position',[700 0 700 1120])
subplot(2,2,1)
plot(SCORE(:,1),SCORE(:,2),'.')
axis tight;
title('PCA1 vs PCA2');
subplot(2,2,2)
plot(LATENT,'x');
axis tight;
title('EigenValues');
subplot(2,2,3:4)
imagesc(hist3([SCORE(1:end,1) SCORE(1:end,2)],[200 200]));
title('PCA1 vs PCA2');

nPC = 5;
num_c = 3;
tt = size(colormap,1);
temp = colormap;
my_colors=temp([1:round(tt/num_c):tt],:);
my_colors = [my_colors; my_colors];
my_style =['.','o','x','+','*','s','d','v','^','<','>','p','h','.','o','x'
    ,'.','o','x','+','*','s','d','v','^','<','>','p','h','.','o','x'];
my_syle = [my_style, my_style, my_style, my_style];

tic
%  [klusters Center]= kmeans(single(allElectrodesXspikes)',num_c);
[klusters Center]= kmeans(SCORE(:,1:nPC),num_c); %% col variables
toc
figure(20)
hold off
set(gcf,'Position',[700 0 1500 1120]);
subplot(2,1,1)
for i =1:num_c
    pp=plot(SCORE(klusters==i,1),SCORE(klusters==i,2),'.');
%     set(pp,'color',my_colors(i,:));
    set(pp,'Marker',my_style(i));
    hold all
    %     i
    %     pause;
end
legend('1', '2', '3' ,'4' ,'5', '6', '7', '8', '9',...
    '10','11','12','13','14','15','16');
title(['Overclustered in PCA1 vs PCA2 (nPC:' num2str(nPC) ')']);

% Combine klusters
% Y = pdist(Center);
Y = pdist(Center,'mahalanobis');
Z = linkage(Y);  % cluster function might be used to autocluster
% dendrogram(Z)
t = .5*(max(Z(:,3)));
subplot(2,1,2);
dendrogram(Z,0,'colorthreshold',t,'orientation','top'); % to show all klusters
%% ADD number of spikes in each cluster
