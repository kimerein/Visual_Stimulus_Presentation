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
colororder =['r.';'g.';'b.';'c.';'m.';'y.';'k.';...
    'ro';'go';'bo';'co';'mo';'yo';'ko';...
    'r.';'g.';'b.';'c.';'m.';'y.';'k.';...
    'ro';'go';'bo';'co';'mo';'yo';'ko';...
    'r.';'g.';'b.';'c.';'m.';'y.';'k.';...
    'ro';'go';'bo';'co';'mo';'yo';'ko';...
    'r.';'g.';'b.';'c.';'m.';'y.';'k.';...
    'ro';'go';'bo';'co';'mo';'yo';'ko';...
    'r.';'g.';'b.';'c.';'m.';'y.';'k.';...
    'ro';'go';'bo';'co';'mo';'yo';'ko'];
Ts = 2.0000e-005;

Z1 = 5;  %Column with Intracellular spike Data 
Z2 = 3;
Y1 = 2;
Y2 = 4;
figure(figure_ind)
subplot(2,1,1)
plot(anal_spikes(:,:,Z1))
subplot(2,1,2)
plot(anal_spikes(:,:,Z2))
graphTitle = '';
%HISTOGRAM
close all;
figure(figure_ind +1)
subplot(2,1,1)
clear Bin_Value; clear Bin_Location;
[Bin_Value Bin_Location]  = hist(anal_spikes(:,:,Z1)',800);
image(Bin_Location,-[1:size(Bin_Value,2)],Bin_Value)
subplot(2,1,2)
clear Bin_Value; clear Bin_Location;
[Bin_Value Bin_Location]  = hist(anal_spikes(:,:,Z2)',800);
image(Bin_Location,-[1:size(Bin_Value,2)],Bin_Value)
text(80,.8,graphTitle,'FontSize',14,'FontWeight','normal','HorizontalAlignment','center' );
ylabel('[mV]')
hold all;


%%% Z1
 [PC, SCORE1, LATENT, TSQUARE] = princomp(anal_spikes(:,:,Z1)');
%  figure(figure_ind +7)
%  plotmatrix(SCORE1(:,1:5))
 %% PLOT SCORE1
 beg_ind = 1;
 end_ind = size(SCORE1,1);
figure(figure_ind +100)
subplot(2,2,1)
plot(SCORE1(beg_ind:end_ind,1),SCORE1(beg_ind:end_ind,2),colororder(1,:))
title('PC1 vs PC2');
hold all;
subplot(2,2,2)
plot(SCORE1(beg_ind:end_ind,1),SCORE1(beg_ind:end_ind,3),colororder(1,:))
title('PC1 vs PC3');
hold all;
subplot(2,2,3)
plot(SCORE1(beg_ind:end_ind,1),SCORE1(beg_ind:end_ind,4),colororder(1,:))
title('PC1 vs PC4');
hold all;
subplot(2,2,4)
plot(SCORE1(beg_ind:end_ind,2),SCORE1(beg_ind:end_ind,3),colororder(1,:))
title('PC2 vs PC3');

%%%%  Z2
[PC, SCORE2, LATENT, TSQUARE] = princomp(anal_spikes(:,:,Z2)');
%  figure(figure_ind +6)
%  plotmatrix(SCORE2(:,1:5))
 %% PLOT SCORE2
 beg_ind = 1;
 end_ind = size(SCORE2,1);
figure(figure_ind +101)
subplot(2,2,1)
plot(SCORE2(beg_ind:end_ind,1),SCORE2(beg_ind:end_ind,2),colororder(1,:))
title('PC1 vs PC2');
hold all;
subplot(2,2,2)
plot(SCORE2(beg_ind:end_ind,1),SCORE2(beg_ind:end_ind,3),colororder(1,:))
title('PC1 vs PC3');
hold all;
subplot(2,2,3)
plot(SCORE2(beg_ind:end_ind,1),SCORE2(beg_ind:end_ind,4),colororder(1,:))
title('PC1 vs PC4');
hold all;
subplot(2,2,4)
plot(SCORE2(beg_ind:end_ind,2),SCORE2(beg_ind:end_ind,3),colororder(1,:))
title('PC2 vs PC3');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NClusters = 2;
% clusters = kmeans(anal_spikes(:,:,Z1)',NClusters);
% %for i=NClusters:-1:1
% graphTitle = 'Intracellular';
% for i=1:NClusters
%     clear Ecluster; clear cluster_ind; clear EcCell;
%     cluster_ind = find(clusters == i); %EXTRACT CLUSTER
%     Ecluster = anal_spikes(:, cluster_ind,Z1);
%     EcCell = anal_spikes(3, cluster_ind,Y1);
%     % %plot of Waveforms from clusters
%     figure(figure_ind+10)
%     plot(Ecluster,colororder(i,:))
%     text(80,.8,graphTitle,'FontSize',14,'FontWeight','normal','HorizontalAlignment','center' );
%     ylabel('[mV]')
%     hold all;
%  %CONTRIBUTION by each cell to group
%     uniq_Cells = unique(EcCell)  ;
%     clear contr_data ;
%     contr_data = zeros(2,size(uniq_Cells,2));
%     graph_text = ''; clear sum_spikes;
%     for ii = 1: size(uniq_Cells,2)
%         sum_spikes(ii) = sum(EcCell == uniq_Cells(ii));%spikes per cell in cluster
%         contr_data(:,ii) = [uniq_Cells(ii) sum_spikes(ii)/size(EcCell,2)]
%         graph_text = sprintf('%s%d\n',graph_text,uniq_Cells(ii));
%     end
%     
%     figure(figure_ind+11)
%     subplot(NClusters,3,i*3-2)
%     plot([1:1:length(Ecluster(:,i))]./Ts,Ecluster)
%     %     text(3e6,.6,graph_text,'FontSize',10,'FontWeight','normal','HorizontalAlignment','left' );    
%     ylabel('[mV]')
%     subplot(NClusters,3,i*3-1)
%     bar(contr_data(2,:));
%     %     set(gca,'XTick',contr_data(1,:))
%     %LIST if cells
%     subplot(NClusters,3,i*3)
%     plot(0,0)
%     graph_text = sprintf('%s\n%d',graph_text,sum(sum_spikes(:)));
%     text(0,0,graph_text,'FontSize',8,'FontWeight','normal','HorizontalAlignment','left' );    
%  
%   %PLOT first 4 PC against each other
%     figure(figure_ind +13)
%     subplot(2,2,1)
%     plot(SCORE1(cluster_ind,1),SCORE1(cluster_ind,2),colororder(i,:))
%     title('PC1 vs PC2');
%     hold all;
%     subplot(2,2,2)
%     plot(SCORE1(cluster_ind,1),SCORE1(cluster_ind,3),colororder(i,:))
%     title('PC1 vs PC3');
%     hold all;
%     subplot(2,2,3)
%     plot(SCORE1(cluster_ind,1),SCORE1(cluster_ind,4),colororder(i,:))
%     title('PC1 vs PC4');
%     hold all;
%     subplot(2,2,4)
%     plot(SCORE1(cluster_ind,2),SCORE1(cluster_ind,3),colororder(i,:))
%     title('PC2 vs PC3');
%     hold all;
%     % figure(figure_ind +9)
%     % for i=1:3
%     % subplot(1,2,i)
%     %     plot(SCORE1(:,i),SCORE1(:,i+1),colororder(i,:))
%     %     title(['PC', num2str(mod(i,3)), 'vs PC', num2str(mod(i+1,3))]);
%     %
%     % hold all;
%     % end
% 
% end
% 
% 
% 
% 
% %EXTRACELLULAR
% 
% figure_ind = figure_ind + 20;
% NClusters = 2;
% clusters = kmeans(anal_spikes(:,:,Z2)',NClusters);
% graphTitle = 'Extracellular';
% alltext = ''; 
% for i=NClusters:-1:1
%     %for i=1:NClusters
%     clear Ecluster; clear cluster_ind; clear EcCell;
%     cluster_ind = find(clusters == i); %EXTRACT CLUSTER
%     Ecluster = anal_spikes(:, cluster_ind,Z2);
%     EcCell = anal_spikes(3, cluster_ind,Y2);  %%%%%%%%%%%%%%%%%
%     uniq_Cells = unique(EcCell)  ;
%     Eplot = find(anal_spikes(3,:,Y2)  == uniq_Cells(1,i));
%     beg_ind = min(Eplot);
%     end_ind = max(Eplot);
%     % %plot of Waveforms from clusters
%     figure(figure_ind+10)
%     plot(Ecluster,colororder(i,:))
%     text(80,2.2,graphTitle,'FontSize',14,'FontWeight','normal','HorizontalAlignment','center' );
%     ylabel('[mV]')
%     hold all;
%     
%  %CONTRIBUTION by each cell to group
%     uniq_Cells = unique(EcCell)  ;
%     clear contr_data ;
%     contr_data = zeros(2,size(uniq_Cells,2));
%     graph_text = ''; clear sum_spikes;
%     for ii = 1: size(uniq_Cells,2)
%         sum_spikes(ii) = sum(EcCell == uniq_Cells(ii));%spikes per cell in cluster
%         contr_data(:,ii) = [uniq_Cells(ii) sum_spikes(ii)/size(EcCell,2)]
%         graph_text = sprintf('%s%d\n',graph_text,uniq_Cells(ii));
%     end
%     figure(figure_ind+11)
%     subplot(NClusters,3,i*3-2)
%     plot([1:1:length(Ecluster(:,i))]./Ts,Ecluster)
%     %     text(3e6,.6,graph_text,'FontSize',10,'FontWeight','normal','HorizontalAlignment','left' );    
%     ylabel('[mV]')
%     subplot(NClusters,3,i*3-1)
%     bar(contr_data(2,:));
%     %     set(gca,'XTick',contr_data(1,:))
%     %LIST if cells
%     subplot(NClusters,3,i*3)
%     plot(0,0)
%     graph_text = sprintf('%s\n%d',graph_text,sum(sum_spikes(:)));
%     text(0,0,graph_text,'FontSize',8,'FontWeight','normal','HorizontalAlignment','left' );    
%  
%    % alltext = sprintf('%s**  %d  **%s',alltext,i,graph_text);
% 
%       %PLOT first 4 PC against each other
%     figure(figure_ind +14)
%     subplot(2,2,1)
%     plot(SCORE2(cluster_ind,1),SCORE2(cluster_ind,2),colororder(i,:))
%     title('PC1 vs PC2');
%     hold all;
%     subplot(2,2,2)
%     plot(SCORE2(cluster_ind,1),SCORE2(cluster_ind,3),colororder(i,:))
%     title('PC1 vs PC3');
%     hold all;
%     subplot(2,2,3)
%     plot(SCORE2(cluster_ind,1),SCORE2(cluster_ind,4),colororder(i,:))
%     title('PC1 vs PC4');
%     hold all;
%     subplot(2,2,4)
%     plot(SCORE2(cluster_ind,2),SCORE2(cluster_ind,3),colororder(i,:))
%     title('PC2 vs PC3');
%     hold all;
%     % figure(figure_ind +9)
%     % for i=1:3
%     % subplot(1,2,i)
%     %     plot(SCORE2(:,i),SCORE2(:,i+1),colororder(i,:))
%     %     title(['PC', num2str(mod(i,3)), 'vs PC', num2str(mod(i+1,3))]);
%     %
%     % hold all;
%     % end
% 
% end

%TODO

%PLOT grouped of spikes in different colors on same plot
%see groups of spikes on PCA space

% are they falling into the same groups?





%{
figure(figure_ind+16)
plotmatrix(SCORE2(cluster1_ind,1:5)) %what does SCORE2 mean?
figure(figure_ind+17)
plotmatrix(SCORE2(cluster2_ind,1:5))

% %
%plot(anal_spikes_Intra(:,:,2));




spike_sort = 1


if spike_sort == 1

figure(4)
plot(my_spikes)
figure(5)
hist(my_ind)

%}