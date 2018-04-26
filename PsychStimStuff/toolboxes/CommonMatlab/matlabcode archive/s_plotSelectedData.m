%% SCRIPT
%% loads file intra  extra cellular data from .mat file
%% plots average of selected files from each file on same axis.
%% DEFAULT file selection: file with largest average negative phase of
%%      extracelluar signal
function [] = plotSelectedData
colororder2 =['r-';'g-';'b-';'c-';'m-';'k-'];
colororder =['r.';'g.';'b.';'c.';'m.';'k.'];

close all;
% readfilename = ['051005_H03.mat';'270905_A02.mat']; % for USER specified
% files
readdirheader = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\mdata\';

readfilename = []; %% for all files
%% all MAT  in directory
temp = strcat(readdirheader,'*.mat');
temp = dir(temp);
for(i=1:size(temp,1))
    if(length(temp(i).name) == 14)
        readfilename = [readfilename; temp(i).name];
    end
end

% savefilename = 'E:\My documents\Academic\Rotation\Scanziani\Data Analysis\ExtraIntra\extracted_spikes.mat'
sselectData = []; selectData = [];selectDataNum = [];
for(i = 1:size(readfilename,1))  %% for all MAT files
    fileIN = strcat(readdirheader,readfilename(i,:));
    [tempselectData tempsselectData tempselectDataNum] = selectMaxNeg(fileIN);
%     [tempselectData tempsselectData tempselectDataNum] = selectAll(fileIN);

    selectDataNum  = [selectDataNum tempselectDataNum];%number of sweeps in mean
    sselectData = [sselectData; tempsselectData];
    selectData = [selectData tempselectData];
end

proData = baselineCorrect(selectData);
proData = normalize(proData);
%% center data on negative peak
offset = center(proData);

Figure_ID = 300;
timedata = [1:1:size(proData,1)];
index = 0; lastCell = [];
%% Plot selected DATA
for (j = 1:size(sselectData,1) )   %%get selected Data from this cells *.mat file
    if (~strcmp(lastCell,sselectData(j,9:11)))%% Color each CEll new color
        lastCell = sselectData(j,9:11);
        index = index +1;
    end
    if (index >  size(colororder,1))
        index = 1;
    end

    figure(Figure_ID)%intra
    %         subplot(2,1,1)
    plot(timedata +60-offset(j),proData(:,j,1),colororder2(index,:));
    line([60,60],[50, -50],'Color','k');
    hold on;
    sSpike= sprintf('%s Sw%d',sselectData(j,:),selectDataNum(1,j));
    title('Max Negative Phase','Interpreter','none');
    text(10,50-j*5,sSpike,'FontSize',8,'Interpreter','none','Color',colororder2(index,1));

    figure(Figure_ID+1)%intra
    %         subplot(2,1,1)
    plot(timedata +60-offset(j),proData(:,j,3),colororder2(index,:));
    line([60,60],[50, -50],'Color','k');
    hold on;
    title('Max Negative Phase','Interpreter','none');
    ylim([-1,2])

end%%END Plot selected DATA

function [selectData sselectData selectDataNum ] = selectMaxNeg(readfilename)
%% Find Largest Negative going extracellular spike
Min =100;
clear indMin;

load(readfilename,...
    'processedFiles','Nspikes_in_File','extracted_spikes', 'EE_coord','-mat');
for iii = 1:size(Nspikes_in_File,2) %% For all data sets stored in file
    start_ind = sum(Nspikes_in_File(1,1:iii-1));
    cData = mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,iii),:),2);
    cMin = min(cData(:,:,3));
    if (cMin < Min)
        Min = cMin;
        indMin = iii;
        Data = cData;
    end
end

temp = sprintf('%s %s',readfilename(end-13:end-4),processedFiles(indMin,12:15));
selectDataNum  = Nspikes_in_File(1,indMin);%number of sweeps in mean
sselectData =  temp;
selectData =  Data;

function [selectData sselectData selectDataNum ] = selectAll(readfilename)
%% Find Largest Negative going extracellular spike
selectDataNum = [];sselectData = []; selectData = [];
load(readfilename,...
    'processedFiles','Nspikes_in_File','extracted_spikes', 'EE_coord','-mat');
for iii = 1:size(Nspikes_in_File,2) %% For all data sets stored in file
    start_ind = sum(Nspikes_in_File(1,1:iii-1));
    cData = mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,iii),:),2);
    temp = sprintf('%s %s',readfilename(end-13:end-4),processedFiles(iii,12:15));
    selectDataNum  = [selectDataNum Nspikes_in_File(1,iii)];%number of sweeps in mean
    sselectData =  [sselectData; temp];
    selectData =  [selectData cData];
end

