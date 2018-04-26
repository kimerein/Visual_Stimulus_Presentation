% % s_kmeanskluster.m
% % use with s_sortybyPCA and s_SpikeSort....
%%
%%  makes num_c kmeans clusters in SCORE 
%% plots and prints to file if breport =1

tic
    %  [klusters Center]= kmeans(single(allElectrodesXspikes)',num_c);
%     [klusters Center]= kmeans(SCORE,num_c,'emptyaction','singleton','distance','correlation');
    [klusters Center]= kmeans(SCORE,num_c,'emptyaction','singleton');
    toc

    if bsortAmp
        kfigID = 21
    else
        kfigID = 20
    end
% end
% Show klusters in 2dim PCA space
figure(kfigID)
hold off
set(gcf,'Position',[700 0 1500 1120]);
        orient tall
if blinkage
    subplot(2,1,1)
end
for i =1:num_c
    pp=plot(SCORE(klusters==i,1),SCORE(klusters==i,2),'.');
    set(pp,'color',my_colors(i,:));
    % set(pp,'Marker',my_style(i));
    hold all
    %     i
    %     pause;
end
legend('1', '2', '3' ,'4' ,'5', '6', '7', '8', '9',...
    '10','11','12','13','14','15','16','17','18','19','20');
if bsortAmp
    title(['Trough vs Peak (' num2str(triggerEE) ')']);
else
    title(['Overclustered in PCA1 vs PCA2 (nPC:' num2str(nPC) ')']);
end
% Combine klusters
if blinkage
    Y = pdist(Center); 
% %     Y = pdist(Center,'seuclidean');
%         Y = pdist(Center,'correlation');  % quite good

    Z = linkage(Y);  % cluster function might be used to autocluster
    % dendrogram(Z)
    t = .5*(max(Z(:,3)));
    subplot(2,1,2);
    dendrogram(Z,0,'colorthreshold',t,'orientation','top'); % to show all klusters
    %% ADD number of spikes in each cluster
end
% if bsave
% end
if breport
      dirtemp = 'REPORT';
      figdesc = ['Kluster_' 'EE' num2str(triggerEE)];
      savefigure(writedirheader,dirtemp,figdesc,'png');  
end