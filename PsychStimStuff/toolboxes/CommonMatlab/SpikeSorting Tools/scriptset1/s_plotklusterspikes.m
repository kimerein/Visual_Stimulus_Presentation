% % s_plotklusterspikes.m
% % use with % % s_kmeanskluster and s_sortybyPCA and s_SpikeSort....

%% plots spike data for each kluster 
NN = 4; % subplots per figure
NP = 5; % plots per cluster;

unit_Nspikes = zeros(1,num_c);
figString = 'KlusterSpikes';
Nwindows = ceil(num_c/NN)
if NN > num_c
    NN = num_c;
end
MAXJ = num_c;
figIND = 200;
str = 'klu';

s_plotLusterSpikes

figure(30)
set(gcf,'Position',[700 0 1500 1120]);
        orient tall
    subplot(2,1,1);

for j=1:size(meanSpike,2)
    if j ==1;
        hold off;
    end
    pp = plot(meanSpike(:,j),'-');
    
%     set(pp,'color',colorO(j,:));
    set(pp,'color',my_colors(j,:));
    hold on;
end
legend('1', '2', '3' ,'4' ,'5', '6', '7', '8', '9',...
    '10','11','12','13','14','15','16','17', '18', '19','20','21','22','23','24','25','26','27', '28', '29','30','31','32','33','34','35','36');
    subplot(2,1,2);
    dendrogram(Z,0,'colorthreshold',t,'orientation','top'); % to show all klusters
