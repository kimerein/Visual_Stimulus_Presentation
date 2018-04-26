%% Plots data from Objective -2 
%%  USED with sO-2_SpikeExtractionAlignPeaks001.m
%% EE and Intra simultaneous recordings
%% EE at different locations in space around the cell

%Plot extracted spikes files in different colors
% selFile_ind = [1:12];
selFile_ind = [1:size(Nspikes_in_File,2)];%for all files
bSamePlot = 0; %plot all in 1 figure or not
b2dplot =1;
b3dplot =1;
bRaw = 1; %plot Raw data

figure_ID = 10;
close all;
colororder2 =['r-';'g-';'b-';'c-';'m-';'y-';'k-'];
colororder =['r.';'g.';'b.';'c.';'m.';'y.';'k.'];

load(savefilename,...
    'processedFiles','Nspikes_in_File','extracted_spikes', '-mat');
index = 1;%for colorindex
timedata = [1:size(extracted_spikes,1)];
for i = 1:size(Nspikes_in_File,2)
    if (~isempty(find(selFile_ind == i)))
        start_ind = sum(Nspikes_in_File(1,1:i-1))
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
                    sSpike= sprintf('%s\n%s\n%d,%d,%d',savefilename(end-20:end),processedFiles(i,6:15),EE_coord(i,2),EE_coord(i,3),EE_coord(i,4));
                    title(sSpike);
                end
                subplot(2,1,2)
                plot(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),3),colororder(index,:));
                hold on;
             
        end
                %plot MEAN data
        figure(100)%intra
        %         subplot(2,1,1)
        plot(mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),1),2),colororder2(index,:));
        hold on;
        %         subplot(2,1,2)
        figure(101)%extra
        plot(mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),3),2),colororder2(index,:));
        hold on;
    
        % PLOT 2D
        if (b2dplot)
            figure(102)%extra
            plot(timedata + EE_coord(i,3)*1,EE_coord(i,2)*1+mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),3),2),colororder2(index,:));
            hold on;
            %Text Box Z= value and Filenumber
            title(savefilename(end-20:end));
            sSpike= sprintf('#%s z=%d',processedFiles(i,6:15),EE_coord(i,4));
            text(EE_coord(i,3)*1,EE_coord(i,2)*1,sSpike,'FontSize',8,'FontWeight','normal','HorizontalAlignment','left' );
        end
        %3D plot
        if(b3dplot)
            Yvalue = EE_coord(i,4)*ones(1,size(extracted_spikes,1));
            figure(103)%extra
            plot3(timedata + EE_coord(i,3)*1,Yvalue,EE_coord(i,2)*1+mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),3),2),colororder2(index,:));
            hold on;
            grid on;
            title(savefilename(end-20:end));
            axis square;
            %         text(EE_coord(i,3)*1,Yvalue,EE_coord(i,2)*1,sSpike,'FontSize',8,'FontWeight','normal','HorizontalAlignment','left' );
        end
        if (bSamePlot ~=1)
            figure_ID = figure_ID +1
        end
        index = index +1;
    end
end