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
    'ro';'go';'bo';'co';'mo';'yo';'ko'];
Ts = 2.0000e-005;
  

% PLOT EXTRA and INTRA SPIKES
figure(figure_ind)
subplot(2,1,1)
plot(extracted_spikes(:,:,5))
subplot(2,1,2)
plot(extracted_spikes(:,:,3))
%plot(extracted_spikes(:,find (extracted_spikes(4,:,4) ==1),3)); 

% add histogram
%% PCA PLOT MATRIX 
[PC, SCORE, LATENT, TSQUARE] = princomp(extracted_spikes(:,:,5)');
figure(figure_ind +7)
plotmatrix(SCORE(:,1:5))
[PC, SCORE, LATENT, TSQUARE] = princomp(extracted_spikes(:,:,3)');
figure(figure_ind +6)
plotmatrix(SCORE(:,1:5))


NClusters = 14;

%% CLUSTER Analysis on INTRA
clusters = kmeans(extracted_spikes(:,:,5)',NClusters);
NPlots = 5;
alltext = '';
for i=1:NPlots
    clear Ecluster; clear cluster_ind; clear EcCell; graph_text = '';
    cluster_ind = find(clusters == i);
    Ecluster = extracted_spikes(:, cluster_ind,5);
    EcCell = extracted_spikes(3, cluster_ind,2);
%CONTRIBUTION by each cell to group
    uniq_Cells = unique(EcCell)  ;
    clear contr_data ;
    contr_data = zeros(2,size(uniq_Cells,2));
    graph_text = ''; clear sum_spikes;
    for ii = 1: size(uniq_Cells,2)
        sum_spikes(ii) = sum(EcCell == uniq_Cells(ii));%spikes per cell in cluster
        contr_data(:,ii) = [uniq_Cells(ii) sum_spikes(ii)/size(EcCell,2)]
        graph_text = sprintf('%s%d\n',graph_text,uniq_Cells(ii));
    end
    figure(figure_ind+15)
    subplot(NPlots,3,i*3-2)
    plot([1:1:length(Ecluster(:,i))]./Ts,Ecluster)
    %     text(3e6,.6,graph_text,'FontSize',10,'FontWeight','normal','HorizontalAlignment','left' );    
    ylabel('[mV]')
    subplot(NPlots,3,i*3-1)
    bar(contr_data(2,:));
    %     set(gca,'XTick',contr_data(1,:))
    %LIST if cells
    subplot(NPlots,3,i*3)
    plot(0,0)
    graph_text = sprintf('%s\n%d',graph_text,sum(sum_spikes(:)));
    text(0,0,graph_text,'FontSize',8,'FontWeight','normal','HorizontalAlignment','left' );    
   alltext = sprintf('%s**  %d  **%s',alltext,i,graph_text);
end
alltext
%figure(figure_ind+11)
%plotmatrix(SCORE(cluster1_ind,1:5))
%figure(figure_ind+12)
%plotmatrix(SCORE(cluster2_ind,1:5))


%NClusters = 2;

%% CLUSTER Extra
clusters = kmeans(extracted_spikes(:,:,3)',NClusters);
for i=1:NPlots
    clear Ecluster; clear cluster_ind; clear EcCell; graph_text = '';
    cluster_ind = find(clusters == i);
    Ecluster = extracted_spikes(:, cluster_ind,3);
    EcCell = extracted_spikes(3, cluster_ind,4);
%CONTRIBUTION by each cell to group
    uniq_Cells = unique(EcCell)  ;
    clear contr_data ;
    contr_data = zeros(2,size(uniq_Cells,2));
    graph_text = ''; clear sum_spikes;
    for ii = 1: size(uniq_Cells,2)
        sum_spikes(ii) = sum(EcCell == uniq_Cells(ii));%spikes per cell in cluster
        contr_data(:,ii) = [uniq_Cells(ii) sum_spikes(ii)/size(EcCell,2)]
        graph_text = sprintf('%s%d\n',graph_text,uniq_Cells(ii));
    end
    figure(figure_ind+11)
    subplot(NPlots,3,i*3-2)
    plot([1:1:length(Ecluster(:,i))]./Ts,Ecluster)
    %     text(3e6,.6,graph_text,'FontSize',10,'FontWeight','normal','HorizontalAlignment','left' );    
    ylabel('[mV]')
    subplot(NPlots,3,i*3-1)
    bar(contr_data(2,:));
    %     set(gca,'XTick',contr_data(1,:))
    %LIST if cells
    subplot(NPlots,3,i*3)
    plot(0,0)
    graph_text = sprintf('%s\n%d',graph_text,sum(sum_spikes(:)));
    text(0,0,graph_text,'FontSize',8,'FontWeight','normal','HorizontalAlignment','left' );    
    alltext = sprintf('%s**  %d  **%s',alltext,i,graph_text);
end
alltext
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