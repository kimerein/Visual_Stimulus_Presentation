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
  
  temp = 1./max(extracted_spikes(:,:,1));
    Scalar =  ones(size(extracted_spikes,1),1)*temp;
    extracted_spikes(:,:,1) = extracted_spikes(:,:,1).*Scalar ;
    temp = -1./min(extracted_spikes(:,:,3));
    Scalar =  ones(size(extracted_spikes,1),1)*temp;
    extracted_spikes(:,:,3) = extracted_spikes(:,:,3).*Scalar ;
    
    
figure(figure_ind)
subplot(2,1,1)
plot(extracted_spikes(:,:,1))
subplot(2,1,2)
plot(extracted_spikes(:,:,3))


%plot(extracted_spikes(:,find (extracted_spikes(4,:,4) ==1),3)); 

% add histogram

[PC, SCORE, LATENT, TSQUARE] = princomp(extracted_spikes(:,:,1)');
figure(figure_ind +5)
plotmatrix(SCORE(:,1:5))
[PC, SCORE, LATENT, TSQUARE] = princomp(extracted_spikes(:,:,3)');
figure(figure_ind +6)
plotmatrix(SCORE(:,1:5))





cluster1 = [];
cluster2 = [];

clusters = kmeans(extracted_spikes(:,:,1)',2);

cluster1_ind = find(clusters == 1);
cluster2_ind = find(clusters == 2);

cluster1 = extracted_spikes(:, cluster1_ind,1); %[wave, cell#]
cluster2 = extracted_spikes(:, cluster2_ind,1);
c1_cell = extracted_spikes(4, cluster1_ind,1+1)
c2_cell = extracted_spikes(4, cluster2_ind,1+1)
%
figure(figure_ind+10)

subplot(2,2,1)
    plot(c1_cell)
subplot(2,2,2)
plot([1:1:length(cluster1(:,1))]./Ts,cluster1)
ylabel('[mV]')
subplot(2,2,3)
    plot(c2_cell)
subplot(2,2,4)
plot([1:1:length(cluster2(:,1))]./Ts,cluster2)
ylabel('[mV]')


%figure(figure_ind+11)
%plotmatrix(SCORE(cluster1_ind,1:5))
%figure(figure_ind+12)
%plotmatrix(SCORE(cluster2_ind,1:5))
cluster1 = [];
cluster2 = [];

clusters = kmeans(extracted_spikes(:,:,3)',3);

cluster1_ind = find(clusters == 1);
cluster2_ind = find(clusters == 2);
cluster3_ind = find(clusters == 3);

cluster1 = extracted_spikes(:, cluster1_ind,3);
cluster2 = extracted_spikes(:, cluster2_ind,3);
cluster3 = extracted_spikes(:, cluster3_ind,3);
c1_cell = extracted_spikes(4, cluster1_ind,3+1)
c2_cell = extracted_spikes(4, cluster2_ind,3+1)
c3_cell = extracted_spikes(4, cluster3_ind,3+1)

%
figure(figure_ind+15)

subplot(3,2,1)
    plot(c1_cell)
subplot(3,2,2)
plot([1:1:length(cluster1(:,1))]./Ts,cluster1)
ylabel('[mV]')
subplot(3,2,3)
    plot(c2_cell)
subplot(3,2,4)
plot([1:1:length(cluster2(:,1))]./Ts,cluster2)
ylabel('[mV]')
subplot(3,2,5)
    plot(c2_cell)
subplot(3,2,6)
plot([1:1:length(cluster3(:,1))]./Ts,cluster3)
ylabel('[mV]')
figure(figure_ind+16)
plotmatrix(SCORE(cluster1_ind,1:5)) %what does SCORE mean?
figure(figure_ind+17)
plotmatrix(SCORE(cluster2_ind,1:5))

% %
%plot(extracted_spikes_Intra(:,:,2));


%{

spike_sort = 1


if spike_sort == 1

figure(4)
plot(my_spikes)
figure(5)
hist(my_ind)

%}