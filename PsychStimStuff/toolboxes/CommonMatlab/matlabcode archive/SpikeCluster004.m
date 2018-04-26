%***********************************************************************
% SpikeCluster.m
%
% Bassam Atallah
% Last change: 8/01/2005
% Reads .atf file, 
% Clusters spikes in 
%                 extracted_spikes(Amplitude,spike#,intra/header/extra/header)
% z = 1 intra, z= 3 extra
% Z = 2 intraheader, z = 4 extraheader
%
%***********************************************************************



% Plot Extracted data
close all;
figure_ind = 50;
colororder =['r.';'g.';'b.';'c.'];
Ts = 2.0000e-005;    
    
figure(figure_ind)
subplot(2,1,1)
plot(extracted_spikes(:,:,5))
subplot(2,1,2)
plot(extracted_spikes(:,:,3))
graphTitle = '';
%HISTOGRAM
close all;
figure(figure_ind +1)
subplot(2,1,1)
clear Bin_Value; clear Bin_Location;
[Bin_Value Bin_Location]  = hist(extracted_spikes(:,:,5)',800);
image(Bin_Location,-[1:size(Bin_Value,2)],Bin_Value)
subplot(2,1,2)
clear Bin_Value; clear Bin_Location;
[Bin_Value Bin_Location]  = hist(extracted_spikes(:,:,3)',800);
image(Bin_Location,-[1:size(Bin_Value,2)],Bin_Value)
text(80,.8,graphTitle,'FontSize',14,'FontWeight','normal','HorizontalAlignment','center' );
ylabel('[mV]')
hold all;



% [PC, SCORE, LATENT, TSQUARE] = princomp(extracted_spikes(:,:,5)');
% figure(figure_ind +7)
% plotmatrix(SCORE(:,1:5))
% [PC, SCORE, LATENT, TSQUARE] = princomp(extracted_spikes(:,:,3)');
% figure(figure_ind +6)
% plotmatrix(SCORE(:,1:5))


NClusters = 2;
Z = 5;
Z1 = 2;
clusters = kmeans(extracted_spikes(:,:,5)',NClusters);
%for i=NClusters:-1:1
graphTitle = 'Intracellular';
for i=1:NClusters
    clear Ecluster; clear cluster_ind; clear EcCell;
cluster_ind = find(clusters == i); %EXTRACT CLUSTER
Ecluster = extracted_spikes(:, cluster_ind,Z); 
EcCell = extracted_spikes(3, cluster_ind,Z1);
% %plot of Waveforms from clusters
figure(figure_ind+10)
plot(Ecluster,colororder(i,:))
text(80,.8,graphTitle,'FontSize',14,'FontWeight','normal','HorizontalAlignment','center' );
ylabel('[mV]')
hold all;

figure(figure_ind+11)
subplot(NClusters,2,i * 2 -1)
plot(EcCell)
subplot(NClusters,2,i * 2)
plot([1:1:length(Ecluster(:,1))]./Ts,Ecluster)
ylabel('[mV]')

 [PC, SCORE, LATENT, TSQUARE] = princomp(extracted_spikes(:,cluster_ind,Z)');
% figure(figure_ind +12)
% plotmatrix(SCORE(:,1:5),colororder(i,:));
% hold all;
%PLOT first 4 PC against each other 
figure(figure_ind +13)
subplot(2,2,1)
plot(SCORE(:,1),SCORE(:,2),colororder(i,:))
title('PC1 vs PC2');
hold all;
subplot(2,2,2)
plot(SCORE(:,1),SCORE(:,3),colororder(i,:))
title('PC1 vs PC3');
hold all;
subplot(2,2,3)
plot(SCORE(:,1),SCORE(:,4),colororder(i,:))
title('PC1 vs PC4');
hold all;
subplot(2,2,4)
plot(SCORE(:,2),SCORE(:,3),colororder(i,:))
title('PC2 vs PC3');
hold all;
% figure(figure_ind +9)
% for i=1:3
% subplot(1,2,i)
%     plot(SCORE(:,i),SCORE(:,i+1),colororder(i,:))
%     title(['PC', num2str(mod(i,3)), 'vs PC', num2str(mod(i+1,3))]);
% 
% hold all;
% end

end




%EXTRACELLULAR

figure_ind = figure_ind + 20;
NClusters = 2;
Z = 3;
Z1 = 4;
clusters = kmeans(extracted_spikes(:,:,Z)',NClusters);
graphTitle = 'Extracellular';
for i=NClusters:-1:1
%for i=1:NClusters
    clear Ecluster; clear cluster_ind; clear EcCell;
cluster_ind = find(clusters == i); %EXTRACT CLUSTER
Ecluster = extracted_spikes(:, cluster_ind,Z); 
EcCell = extracted_spikes(3, cluster_ind,Z1);  %%%%%%%%%%%%%%%%% 
% %plot of Waveforms from clusters
figure(figure_ind+10)
plot(Ecluster,colororder(i,:))
text(80,2.2,graphTitle,'FontSize',14,'FontWeight','normal','HorizontalAlignment','center' );
ylabel('[mV]')
hold all;
%% Histogram plot

figure(figure_ind+11)
subplot(NClusters,2,i * 2 -1)
plot(EcCell)
subplot(NClusters,2,i * 2)
plot([1:1:length(Ecluster(:,1))]./Ts,Ecluster)
ylabel('[mV]')

 [PC, SCORE, LATENT, TSQUARE] = princomp(extracted_spikes(:,cluster_ind,Z)');
% figure(figure_ind +12)
% plotmatrix(SCORE(:,1:5),colororder(i,:));
% hold all;
%PLOT first 4 PC against each other 
figure(figure_ind +14)
subplot(2,2,1)
plot(SCORE(:,1),SCORE(:,2),colororder(i,:))
title('PC1 vs PC2');
hold all;
subplot(2,2,2)
plot(SCORE(:,1),SCORE(:,3),colororder(i,:))
title('PC1 vs PC3');
hold all;
subplot(2,2,3)
plot(SCORE(:,1),SCORE(:,4),colororder(i,:))
title('PC1 vs PC4');
hold all;
subplot(2,2,4)
plot(SCORE(:,2),SCORE(:,3),colororder(i,:))
title('PC2 vs PC3');
hold all;
% figure(figure_ind +9)
% for i=1:3
% subplot(1,2,i)
%     plot(SCORE(:,i),SCORE(:,i+1),colororder(i,:))
%     title(['PC', num2str(mod(i,3)), 'vs PC', num2str(mod(i+1,3))]);
% 
% hold all;
% end

end

%TODO

%PLOT grouped of spikes in different colors on same plot
%see groups of spikes on PCA space

% are they falling into the same groups?





%{
figure(figure_ind+16)
plotmatrix(SCORE(cluster1_ind,1:5)) %what does SCORE mean?
figure(figure_ind+17)
plotmatrix(SCORE(cluster2_ind,1:5))

% %
%plot(extracted_spikes_Intra(:,:,2));




spike_sort = 1


if spike_sort == 1

figure(4)
plot(my_spikes)
figure(5)
hist(my_ind)

%}