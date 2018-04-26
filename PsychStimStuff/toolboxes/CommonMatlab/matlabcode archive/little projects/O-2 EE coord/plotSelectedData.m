%% SCRIPT
%% Objective -2
%% loads file intra  extra cellular data from .mat file
%% plots average of selected files from each file on same axis.
%% DEFAULT file selection: file with largest average negative phase of
%%      extracelluar signal
%% Extract data with SpikeExtractionAlignPeaks002.m

function [] = plotSelectedData
colororder2 =['r-';'g-';'b-';'c-';'m-';'k-'];
colororder =['r.';'g.';'b.';'c.';'m.';'k.'];

close all;
% readfilename = ['051005_H03.mat';'270905_A02.mat']; % for USER specified
% files
readdirheader = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Data Analysis\mdata\';
% readdirheader = 'E:\My documents\Scanziani Lab\Data Analysis\mdata\';

readfilename = []; %% for all files
%% all MAT  in directory
temp = strcat(readdirheader,'*typef*.mat');
filedata = struct2cell(dir(temp));

if isempty(filedata)
    stemp = sprintf('No Files found at: %s\nCheck path',temp)
    error(stemp);
end
% savefilename = 'E:\My documents\Academic\Rotation\Scanziani\Data Analysis\ExtraIntra\extracted_spikes.mat'
sselectData = cell(0); selectData = [];selectDataNum = [];selectTs=[];
for(i = 1:size(filedata,2))  %% for all MAT files
    readfilename =   filedata{1,i} ;
    fileIN = strcat(readdirheader,readfilename);
        [tempselectData tempsselectData tempselectDataNum Ts tmeanmisaligned] = selectMaxNeg(fileIN);
%       DOESN't work propery need to deal with multiple files in [tempselectData tempsselectData tempselectDataNum Ts] = selectAll(fileIN);


    if isnan(Ts)
        warning('Ts does not exist in .mat file default value 1/0.02e-3 used')
        Ts = 1/0.02e-3;
    end
    selectDataNum  = [selectDataNum tempselectDataNum];%number of sweeps in mean
    sselectData{i,1} = tempsselectData;
    selectData{i,1} =  tempselectData;
    selectTs = [selectTs Ts] ;
%     meanmisaligned(i) = tmeanmisaligned;
end

proData = baselineCorrect(selectData);
proData = normalize(proData);
%% center data on negative peak
 offset = center(proData);
%% align Intra
offset2 = centerIntra(proData);

Figure_ID = 300;
timedata = [1:1:size(proData,1)];
index = 0; lastCell = [];
%% Plot selected DATA
tempTs =  max(selectTs(:));  %need because there may be different Ts for different files
for (j = 1:size(sselectData,1) )   %%get selected Data from this cells *.mat file
    timedata = [1:size(proData{j,1},1)];
    if (~strcmp(lastCell,sselectData{j,1}(strfind(sselectData{j,1}(1,:),'.mat')-3:strfind(sselectData{j,1}(1,:),'.mat')-1)))%% Color each CEll new color
        lastCell = sselectData{j,1}(9:11);
        index = index +1;
    end
    if (index >  size(colororder,1))
        index = 1;
    end

    figure(Figure_ID)%intra
    %         subplot(2,1,1)
    plot((timedata -offset2(j))/selectTs(j),proData{j,1}(:,1,1),colororder2(index,:));
    %     line([60,60],[50, -50],'Color','k');
    hold on;
    sSpike= sprintf('%s Sw%d',sselectData{j,1},selectDataNum(1,j));
    title('Max Negative Phase','Interpreter','none');
        text(45/tempTs,-j*5,sSpike,'FontSize',8,'Interpreter','none','Color',colororder2(index,1));

    figure(Figure_ID+1)%extra
    %         subplot(2,1,1)
    plot((timedata -offset(j))/selectTs(j),proData{j,1}(:,1,3),colororder2(index,:));
    %     line([60,60],[50, -50],'Color','k');
    hold on;
         text(45/tempTs,0-j*.07,sSpike,'FontSize',8,'Interpreter','none','Color',colororder2(index,1));

    title('Max Negative Phase','Interpreter','none');
    %     ylim([-1,2])

    pause;
end%%END Plot selected DATA

function [selectData sselectData selectDataNum Ts meanmisaligned] = selectMaxNeg(readfilename)
%% Find Largest Negative going extracellular spike
Min =100;
clear indMin;

load(readfilename,...
    'processedFiles','Ts','Nspikes_in_File','extracted_spikes', 'misalign_ind', 'EE_coord','-mat');
if ~exist('Ts')
    Ts = 1/20e-6;
end
for iii = 1:size(Nspikes_in_File,2) %% For all data sets stored in file
    start_ind = sum(Nspikes_in_File(1,1:iii-1));
    cData = mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,iii),:),2);
    if ~isempty(cData)
        cMin = min(cData(:,:,3));
        if (cMin < Min)
            Min = cMin;
            indMin = iii;
%             meanmisaligned =  mean(misalign_ind(start_ind+1:start_ind+Nspikes_in_File(1,iii)));
        end
    end
end

% EXCLUDE spikes that that differ from the MODE
start_ind = sum(Nspikes_in_File(1,1:indMin-1));
[spike_ind xcl_ind] = xcludeSpikes(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,indMin),1),Ts);
Data = mean(extracted_spikes(:,spike_ind+start_ind,:),2);

sfilename = readfilename(max(strfind(readfilename,'\'))+1:end)
% deal with fact that some processedFiles are cells not strings
if iscell(processedFiles)
    temp = processedFiles{indMin,1}(end-8:end-4);
else
    temp = processedFiles(indMin,end-8:end-4);
end
stemp =sqrt(EE_coord(indMin,2)^2+ EE_coord(indMin,3)^2 + EE_coord(indMin,4)^2);%distance of EE from cell
temp = sprintf('%s %s %3.1fum',sfilename,temp,stemp );
selectDataNum  = size(spike_ind);%number of sweeps in mean
sselectData =  temp;
selectData =  Data;
if ~exist('Ts')
    Ts = NaN;
end
meanmisaligned =[];
% 
