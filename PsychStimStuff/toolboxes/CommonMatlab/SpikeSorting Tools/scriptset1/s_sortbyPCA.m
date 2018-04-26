% PCA- kmeans sort and plot
%s_sortbyPCA
% CALL from with s_SpikeSort....
    % n_units = num_c;

    % PCA
    if bsortAmp
     [PC, SCORE, LATENT, TSQUARE] = princomp(single(parameters));
    else
     [PC, SCORE, LATENT, TSQUARE] = princomp(single(allElectrodesXspikes'));
    end
    % end
    % Scatter plot first two PCs
    figure(10+bsortAmp)
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

    if breport
        dirtemp = 'REPORT';
        figdesc = ['PCA12_' 'EE' num2str(triggerEE)];
        savefigure(writedirheader,dirtemp,figdesc,'pdf');  
    end
    clear PC; clear LATENT; clear TSQUARE;

