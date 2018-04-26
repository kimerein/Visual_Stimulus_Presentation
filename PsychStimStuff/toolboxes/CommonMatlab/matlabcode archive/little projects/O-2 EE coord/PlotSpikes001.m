%% Plots data from Objective -2 
%%  USED with sO-2_SpikeExtractionAlignPeaks001.m
%% EE and Intra simultaneous recordings
%% EE at different locations in space around the cell

%Plot extracted spikes files in different colors
% selFile_ind = [1:12];
selFile_ind = [1:size(Nspikes_in_File,2)];%for all files
bSamePlot = 0; %plot all in 1 figure or not
b2dplot =0;
b3dplot =1;
bRaw = 1; %plot Raw data

figure_ID_init = 10;
figure_ID = figure_ID_init;
close all;
colororder2 =['r-';'g-';'b-';'c-';'m-';'y-';'k-'];
backcolororder2 =['g-';'b-';'c-';'m-';'y-';'k-';'r-'];
colororder =['r.';'g.';'b.';'c.';'m.';'y.';'k.'];
backcolororder =['g.';'b.';'c.';'m.';'y.';'k.';'r.'];

savefileheader = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Data Analysis\mdata\'; 
stemp = strcat(savefileheader,savefilename);
load(stemp,...
    'Ts','processedFiles','Nspikes_in_File','extracted_spikes', 'misalign_ind','EE_coord','-mat');
if (~iscell(processedFiles)) %% fix some older files where processFiles was not a cell
    temp = cell(size(processedFiles,1),1);
    for i=1:size(processedFiles,1)
        temp{i,1} = processedFiles(i,:);
    end
    processedFiles = temp;
%   save(savefilename,...
%     'Ts','processedFiles','Nspikes_in_File','extracted_spikes', 'EE_coord','-mat');
end
sstemp = '';
index = 1;%for colorindex
timedata = [1:size(extracted_spikes,1)]/Ts;
for i = 1:size(Nspikes_in_File,2)
    if (~isempty(find(selFile_ind == i)))
        start_ind = sum(Nspikes_in_File(1,1:i-1))
        if (index >  size(colororder,1))
            index = 1;
        end
        
                % exclude spikes that occur too soon
        [spike_ind xcl_ind] = xcludeSpikes(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),1),Ts);
        %         plot RAW data


        if (bRaw)
                figure(figure_ID)%intra
                subplot(2,1,1)
                plot(timedata,extracted_spikes(:,start_ind+spike_ind,1),colororder(index,:));
                if ~isempty(xcl_ind)
                    hold on;
                  plot(timedata,extracted_spikes(:,start_ind + xcl_ind,1),'-k');
                end
                hold on;
                if (bSamePlot~=1)
                    sstemp = sprintf('%s %s,',sstemp,processedFiles{i,1}(6:15));
                    sstemp2 = savefilename(max(strfind(savefilename,'\'))+1:end);
                    sSpike= sprintf('%s\n%s\n%d,%d,%d',sstemp2,processedFiles{i,1}(6:15),EE_coord(i,2),EE_coord(i,3),EE_coord(i,4));
                    title(sSpike,'Interpreter','none');
                end
                subplot(2,1,2)
                plot(timedata,extracted_spikes(:,start_ind+spike_ind,3),colororder(index,:));
                if ~isempty(xcl_ind)
                    hold on;
                    plot(timedata,extracted_spikes(:,start_ind+xcl_ind,3),'-k');
                end
                hold on;
             
        end

        meanIntra = mean(extracted_spikes(:,start_ind+spike_ind,1),2);
        meanExtra = mean(extracted_spikes(:,start_ind+spike_ind,3),2);
  
                %plot MEAN data
        figure(100)%intra
        %         subplot(2,1,1)
        plot(timedata ,meanIntra,colororder2(index,:));
        sSpike= sprintf('%s\n%s',sstemp2,sstemp); title(sSpike,'Interpreter', 'none');
        hold on;
        xlim([0 max(timedata)]);
        %         subplot(2,1,2)
        figure(101)%extra
        
        plot(timedata,meanExtra,colororder2(index,:));
        sSpike= sprintf('%s\n%s',sstemp2,sstemp); title(sSpike,'Interpreter', 'none');
        hold on;
        xlim([0 max(timedata)]);

    
        % PLOT 2D
        if (b2dplot)
            figure(102)%extra
            plot(timedata + EE_coord(i,3)*1*1/Ts,EE_coord(i,2)*1+meanExtra,colororder2(index,:));
            hold on;
            %Text Box Z= value and Filenumber
            title(sstemp2,'Interpreter','none');
            sSpike= sprintf('#%s z=%d',processedFiles{i,1}(6:15),EE_coord(i,4));
            text(EE_coord(i,3)*1,EE_coord(i,2)*1,sSpike,'FontSize',8,'FontWeight','normal','HorizontalAlignment','left','Interpreter','none');
        end
        %3D plot
        if(b3dplot)
            Yvalue = EE_coord(i,4)*ones(1,size(extracted_spikes,1));
            figure(103)%extra
            plot3(timedata + EE_coord(i,3)*1,Yvalue,EE_coord(i,2)*1+meanExtra,colororder2(index,:));
            hold on;
            grid on;
            title(savefilename,'Interpreter','none');
            axis square;
            %         text(EE_coord(i,3)*1,Yvalue,EE_coord(i,2)*1,sSpike,'FontSize',8,'FontWeight','normal','HorizontalAlignment','left' );
        end
        if (bSamePlot ~=1)
            figure_ID = figure_ID +1
        end
        index = index +1;
    end
end

if false
    %% PRINT DATA TO FILE
    savefileheader = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Data Analysis\Figures\';
    if strfind(savefilename,'\')
        filename = savefilename(max(strfind(savefilename,'\'))+1:end)
    else
        filename = savefilename;
    end
    savepicfilename = strcat(savefileheader,filename);

    for i=0:(index-2)
        figure(figure_ID_init + i)
        temp = sprintf('%sRawData%d.tif',savepicfilename(1:end-4),i+1);
        print('-dtiffn',temp)
    end

    figure(100)
    temp = sprintf('%sMeanIntraSpike.tif',savepicfilename(1:end-4));
    print('-dtiffn',temp)

    figure(101)
    temp = sprintf('%sMeanExtraSpike.tif',savepicfilename(1:end-4));
    print('-dtiffn',temp)


end
