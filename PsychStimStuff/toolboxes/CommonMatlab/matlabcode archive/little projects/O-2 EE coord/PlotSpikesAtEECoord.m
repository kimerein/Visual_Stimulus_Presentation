%Plot extracted spikes files in different colors
selFile_ind = [1:12];
bSamePlot = 1; %plot all in 1 figure or not
figure_ID = 10;
close all;
colororder2 =['r-';'g-';'b-';'c-';'m-';'y-';'k-'];
colororder =['r.';'g.';'b.';'c.';'m.';'y.';'k.'];

load(savefilename,...
    'EE_coord','processedFiles','Nspikes_in_File','extracted_spikes', '-mat');
index = 1;%for colorindex
timedata = [1:size(extracted_spikes,1)];
for i = 1:size(Nspikes_in_File,2)
    if (~isempty(find(selFile_ind == i)))
        start_ind = sum(Nspikes_in_File(1,1:i-1))
        if (index >  size(colororder,1))
            index = 1;
        end
%         %plot MEAN data
        figure(figure_ID)%intra
%         subplot(2,1,1)
        plot(mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),1),2),colororder2(index,:));
        hold on;
%         subplot(2,1,2)
        figure(figure_ID+1)%extra
        plot(timedata + EE_coord(i,3)*1,EE_coord(i,2)*1+mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),3),2),colororder2(index,:));
        hold on;
        %Text Box Z= value and Filenumber
        title(savefilename(end-20:end));
        sSpike= sprintf('#%s z=%d',processedFiles(i,6:15),EE_coord(i,z));
        text(EE_coord(i,3)*1,EE_coord(i,2)*1,sSpike,'FontSize',8,'FontWeight','normal','HorizontalAlignment','left' );
 %3D plot
 Yvalue = EE_coord(i,4)*ones(1,size(extracted_spikes,1));
        figure(figure_ID+2)%extra
        plot3(timedata + EE_coord(i,3)*1,Yvalue,EE_coord(i,2)*1+mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),3),2),colororder2(index,:));
        hold on;
        grid on;
        axis square;
%         text(EE_coord(i,3)*1,Yvalue,EE_coord(i,2)*1,sSpike,'FontSize',8,'FontWeight','normal','HorizontalAlignment','left' );
 

        if (bSamePlot ~=1)
            figure_ID = figure_ID +1
        end
        index = index +1;
    end
end