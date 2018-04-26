%% SCRIPT
%% Plots abf files with specified header and filenumber
%% xtract Spikes & saves at specified location

%% USER must enter:
%%    readfilenumber, EE_coord
% readfilenumber = []
% EE_coord = [] %X,Y,Z coordinates of extracellular electrode
readfilenumber  = 2;
close all;
readfileheader = ['2005_10_21_'];
savefilename = ['102105_O001_E001.mat'];
savedirheader = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\mdata\';
% savefilename = 'E:\My documents\Academic\Rotation\Scanziani\Data Analysis\ExtraIntra\extracted_spikes.mat'
savefilename = strcat(savedirheader,savefilename);

[processedFiles,Nspikes_in_File,extracted_spikes] = xSpikes(readfileheader,readfilenumber);

save(savefilename,...
    'processedFiles','Nspikes_in_File','extracted_spikes', 'EE_coord','-mat');

%%%%%

% savefilename = ['280905_C03.mat'];
% savedirheader = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\mdata\';
% savefilename = 'E:\My documents\Academic\Rotation\Scanziani\Data Analysis\ExtraIntra\extracted_spikes.mat'
% savefilename = strcat(savedirheader,savefilename);

% close all;
ee3dplot(savefilename,[]);

% plotabf('2005_09_30_',3);
