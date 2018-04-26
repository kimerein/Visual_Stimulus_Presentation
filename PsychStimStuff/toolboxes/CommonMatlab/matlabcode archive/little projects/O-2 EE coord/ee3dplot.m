function [] = ee3dplot(readfilename,selFile_ind)
% Function plots Extracellular spikes extracted by xSpike.m
% 2d and 3d plots are made using the location of the EE in EE_coord

% selFile_ind contains the indices of the abf files that will be plotted
% selFile_ind = [] plots data from all abf files in readfilename

%toggle plots
bSamePlot = 0; %plot all in 1 figure or not
b2dplot =0;
b3dplot =1;
bRaw = 1; %plot spikes aligned to Intra Peak before averaging
%%%
colororder2 =['r-';'g-';'b-';'c-';'m-';'y-';'k-'];
colororder =['r.';'g.';'b.';'c.';'m.';'y.';'k.'];
figure_ID = 10;
load(readfilename,...
    'processedFiles','Nspikes_in_File','extracted_spikes','EE_coord','-mat');
if isempty(selFile_ind)
    selFile_ind = [1:size(Nspikes_in_File,2)];%for all files
end

index = 1;%for colorindex
timedata = [1:size(extracted_spikes,1)];
for i = 1:size(Nspikes_in_File,2)
    if (~isempty(find(selFile_ind == i)))
        start_ind = sum(Nspikes_in_File(1,1:i-1));
        if (index >  size(colororder,1))
            index = 1;
        end
        %         plot RAW data
        if (bRaw)
                figure(figure_ID)%intra
                subplot(2,1,1)
                plot(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),1),colororder(index,:));
                hold on;
                if (bSamePlot~=1)
                    sSpike= sprintf('%s\n%s\n%d,%d,%d',readfilename(end-20:end),processedFiles(i,6:15),EE_coord(i,2),EE_coord(i,3),EE_coord(i,4));
                    title(sSpike,'Interpreter','none');
                end
                subplot(2,1,2)
                plot(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),3),colororder(index,:));
                hold on;
             
        end
    %plot MEAN data
        figure(99)%intra
        %         subplot(2,1,1)
        plot(mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),1),2),colororder2(index,:));
        line([40,40],[50, -50],'Color','k');
        hold on;
        sSpike= sprintf('%s',readfilename(end-20:end));
        title(sSpike,'Interpreter','none');
        sSpike = sprintf('%s',processedFiles(i,12:15));
        text(10,50-i*5,sSpike,'FontSize',8,'Interpreter','none','Color',colororder2(index,1));
        ylim([-70,50])

        %         subplot(2,1,2)
        figure(100)%extra
        plot(mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),3),2)+i*10,colororder2(index,:));
        line([40,40],[10*i, -10],'Color','k');
        hold on;
        sSpike= sprintf('%s',readfilename(end-20:end));
        title(sSpike,'Interpreter','none');
 %         ylim([-20,20])
        figure(101)%extra
        plot(mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),3),2),colororder2(index,:));
        line([40,40],[20, -50],'Color','k');
        hold on;
        sSpike= sprintf('%s',readfilename(end-20:end));
        title(sSpike,'Interpreter','none');
        ylim([-70,30])
        % PLOT 2D
        if (b2dplot)
            figure(102)%extra
            plot(timedata + EE_coord(i,3)*1,EE_coord(i,2)*1+mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),3),2),colororder2(index,:));
            hold on;
            %Text Box Z= value and Filenumber
            title(readfilename(end-20:end));
            sSpike= sprintf('#%s z=%d',processedFiles(i,6:15),EE_coord(i,4));
            text(EE_coord(i,3)*1,EE_coord(i,2)*1,sSpike,'FontSize',8,'FontWeight','normal','HorizontalAlignment','left', 'Interpreter', 'none','Color',colororder2(index,1) );
        end
        %3D plot
        if(b3dplot)
            Yvalue = EE_coord(i,4)*ones(1,size(extracted_spikes,1));
            figure(103)%extra
            plot3(timedata + EE_coord(i,2)*1,Yvalue,EE_coord(i,3)*1+mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),3),2),colororder2(index,:));
            hold on;
            grid on;
            title(readfilename(end-20:end), 'Interpreter', 'none');
            axis square;
            %         text(EE_coord(i,3)*1,Yvalue,EE_coord(i,2)*1,sSpike,'FontSize',8,'FontWeight','normal','HorizontalAlignment','left' );
        end
        if (bSamePlot ~=1)
            figure_ID = figure_ID +1;
        end
        index = index +1;
    end
end