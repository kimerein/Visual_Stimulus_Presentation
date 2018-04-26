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
    %Scale Amplitude to 1
  %{
  temp = 1./max(extracted_spikes(:,:,1));
    Scalar =  ones(size(extracted_spikes,1),1)*temp;
    extracted_spikes(:,:,1) = extracted_spikes(:,:,1).*Scalar ;
    %}
    temp = -1./min(extracted_spikes(:,:,3));
    Scalar =  ones(size(extracted_spikes,1),1)*temp;
    extracted_spikes(:,:,3) = extracted_spikes(:,:,3).*Scalar ;
    
    
figure(figure_ind)
subplot(2,1,1)
plot(extracted_spikes(:,:,5))
subplot(2,1,2)
plot(extracted_spikes(:,:,3))


%plot(extracted_spikes(:,find (extracted_spikes(4,:,4) ==1),3)); 

% add histogram

[PC, SCORE, LATENT, TSQUARE] = princomp(extracted_spikes(:,:,5)');
figure(figure_ind +7)
plotmatrix(SCORE(:,1:5))
[PC, SCORE, LATENT, TSQUARE] = princomp(extracted_spikes(:,:,3)');
figure(figure_ind +6)
plotmatrix(SCORE(:,1:5))


NClusters = 2;

clusters = kmeans(extracted_spikes(:,:,5)',NClusters);

for i=1:NClusters
    clear Ecluster; clear cluster_ind; clear EcCell;
cluster_ind = find(clusters == i);
Ecluster = extracted_spikes(:, cluster_ind,5); 
EcCell = extracted_spikes(3, cluster_ind,2);
%
figure(figure_ind+10)

subplot(NClusters,2,i * 2 -1)
plot(EcCell)
subplot(NClusters,2,i * 2)
plot([1:1:length(Ecluster(:,1))]./Ts,Ecluster)
ylabel('[mV]')

end

%figure(figure_ind+11)
%plotmatrix(SCORE(cluster1_ind,1:5))
%figure(figure_ind+12)
%plotmatrix(SCORE(cluster2_ind,1:5))


%NClusters = 2;

clusters = kmeans(extracted_spikes(:,:,3)',NClusters);

for i=1:NClusters
    clear Ecluster; clear cluster_ind; clear EcCell;
cluster_ind = find(clusters == i);
Ecluster = extracted_spikes(:, cluster_ind,3); 
EcCell = extracted_spikes(3, cluster_ind,3+1);
%
figure(figure_ind+15)

subplot(NClusters,2,i * 2 -1)
plot(EcCell)
subplot(NClusters,2,i * 2)
plot([1:1:length(Ecluster(:,1))]./Ts,Ecluster)
ylabel('[mV]')

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