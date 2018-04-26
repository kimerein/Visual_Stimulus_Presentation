% extracted_spikes_OUT = extracted_spikes(1:110,:,:);
% clear extracted_spikes;
% extracted_spikes = extracted_spikes_OUT;
%
close all;
clear str_legend ;
figure_ind = 100;
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
Z1 = 5;
Y1 = 2;
[PC1, SCORE1, LATENT1, TSQUARE1] = princomp(anal_spikes(:,:,Z1)');
Z2 = 3;
Y2 = 4; %% not used assumed to be the same as Y1
[PC2, SCORE2, LATENT2, TSQUARE2] = princomp(anal_spikes(:,:,Z2)');

graphTitle = 'cell by cell';
EcCell = anal_spikes(3, :,Y1);
uniq_Cells = unique(EcCell)  ;
for i=1 : size(uniq_Cells,2)
    Eplot = find(anal_spikes(3,:,Y1)  == uniq_Cells(1,i));
    %% Std of first 30points (call it noise)
    %% (Max - Min) (signal)

    noise(i) = mean(std(anal_spikes(1:22,Eplot,Z1)));
    signal(i)= mean(max(anal_spikes(:,Eplot,Z1)) - min(anal_spikes(:,Eplot,Z1)));

    figure(figure_ind+1)
    subplot(2,1,1)
    plot_handle = plot(anal_spikes(:,Eplot,Z1),colororder(i,:))
    %text(80,60,graphTitle,'FontSize',14,'FontWeight','normal','HorizontalAl
    %ignment','center' );
    ylabel('[mV]')
    for ii = 1 : i
        temp = sprintf('%2d * %s', ii,num2str(uniq_Cells(1,ii)))
        if length(temp) < 15
            temp = [temp blanks(15-length(temp))];
        end
        str_legend(ii,:) = temp;
    end
    legend(plot_handle,str_legend,'Location','NorthWest')
    hold all;
    subplot(2,1,2)
    plot(anal_spikes(:,Eplot,Z2),colororder(i,:))
    %  text(80,2.2,graphTitle,'FontSize',14,'FontWeight','normal','HorizontalAlignment','center' );
    ylabel('[mV]')
    hold all;

    %FIND Ind of  cell_i in SCORE

    beg_ind = min(Eplot);
    end_ind = max(Eplot);
    % PLOT SCORE1
    figure(figure_ind +10)
    subplot(2,2,1)
    plot(SCORE1(beg_ind:end_ind,1),SCORE1(beg_ind:end_ind,2),colororder(i,:))
    title('PC1 vs PC2');
    hold all;
    subplot(2,2,2)
    plot(SCORE1(beg_ind:end_ind,1),SCORE1(beg_ind:end_ind,3),colororder(i,:))
    title('PC1 vs PC3');
    hold all;
    subplot(2,2,3)
    plot(SCORE1(beg_ind:end_ind,1),SCORE1(beg_ind:end_ind,4),colororder(i,:))
    title('PC1 vs PC4');
    hold all;
    subplot(2,2,4)
    plot(SCORE1(beg_ind:end_ind,2),SCORE1(beg_ind:end_ind,3),colororder(i,:))
    title('PC2 vs PC3');
    hold all;

    %% PLOT SCORE2
    figure(figure_ind +11)
    subplot(2,2,1)
    plot(SCORE2(beg_ind:end_ind,1),SCORE2(beg_ind:end_ind,2),colororder(i,:))
    title('PC1 vs PC2');
    hold all;
    subplot(2,2,2)
    plot(SCORE2(beg_ind:end_ind,1),SCORE2(beg_ind:end_ind,3),colororder(i,:))
    title('PC1 vs PC3');
    hold all;
    subplot(2,2,3)
    plot(SCORE2(beg_ind:end_ind,1),SCORE2(beg_ind:end_ind,4),colororder(i,:))
    title('PC1 vs PC4');
    hold all;
    subplot(2,2,4)
    plot(SCORE2(beg_ind:end_ind,2),SCORE2(beg_ind:end_ind,3),colororder(i,:))
    title('PC2 vs PC3');
    hold all;

    pause;
end
figure(figure_ind)
plot(1:size(signal,2),signal./noise,colororder(1,:));
