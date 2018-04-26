%***********************************************************************
%  USER MUST ENTER EE_coord by hand
%
%  
% Based on SpikeExtraction008
%
% Bassam Atallah
% Last change: 9/29/2005
% Reads .abf file, 
% xSpikes placed in
%                 extracted_spikes(Amplitude,spike#,intra/header/extra/header)
% z = 1 intra, z= 3 extra
% Z = 2 intraheader, z = 4 extraheader
%
%***********************************************************************
% Initialize workspace

clc;
% warning  off all;

%INPUT File parameters
% dirheader =  'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Patch Data\'; 
% dirheader =  'e:\My Documents\Scanziani Lab\Patch Data\'; %HOME
% dirheader =  '\\132.239.158.164\Patch Data\';
dirheader =  'L:\Patch Data\'; %HOME
% dirheader = 'Z:\Patch Data\Modified\';
% dirheader =  '\\132.239.158.164\Patch Data\';
%  dirheader =  '\\Lindseysrig\My Documents\Patch Data\';
% savefileheader = 'E:\My documents\Scanziani Lab\Data Analysis\mdata\';
% savefileheader = 'Z:\Data Analysis\mdata\';
savefileheader = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Data Analysis\mdata\';

% USER ENTER 
[fileheader filenumbers sSliceCellID] = readaxonDialog(1);
Nfiles_per_header = size(filenumbers,2);
% fileheader = ['2005_09_28_'];
% filenumbers = [10:11 13:22]'; %file index for day
% filenumbers = [56];
% Nfiles_per_header = [5];


% CREATe Savefile
temp = strfind(fileheader,'_');
stemp = strcat(fileheader(temp(1)+1:temp(1)+2),fileheader(temp(2)+1:temp(2)+2),fileheader(temp(1)-2:temp(1)-1));
savefilename = sprintf('O-2_%s_%s.mat',stemp,sSliceCellID);

%USE ONLY first time
% processedFiles = char(zeros(1,size(filename,2)));
processedFiles = cell( size(Nfiles_per_header,2),1);
putative_type = [];
cell_number = [];
ID_SpikeLast = 0;
total_extractedspikes = 0;
extracted_spikes = [];
Nspikes_in_File = [];

if(sum(Nfiles_per_header) ~= size(filenumbers,2) & sum(Nfiles_per_header) ~= size(filenumbers,1) & sum(Nfiles_per_header) ~= size(CELL_NUM) & sum(Nfiles_per_header) ~= size(PUT_TYP))
    Warning = 'CHECK INPUT File parameters   Nfiles_per_header ~= size(filenumbers,2)'
end
i = 0;
nExtra_Chn = 0;

% load(strcat(savefileheader,savefilename),'EE_coord', '-mat')
%//////////////////////////
index_file = 0;
for (j = 1: size(fileheader,1))
   filenumber_ind = i+1;
    for (i =filenumber_ind:filenumber_ind + Nfiles_per_header(j,1) -1)  % count through filenumbers for this fileheader
        %Set Filename
        filename = createaxonfilename(fileheader(j,:),filenumbers(i))
        file = strcat(dirheader,filename) ;
        Ts = 1/0.02e-3 ;%Sample rate 50kHz
        my_data_zind = 1; % DON"T KNOW why this variable exists
  
        %CHECK if same file is already loaded (saves time)
clear        lastfileloaded;
        b_loadnewdata =1
        if(exist('lastfileloaded'))
            if (strcmp(file,lastfileloaded))
                b_loadnewdata = 0;
            end
        end
b_loadnewdata = 1;
        if b_loadnewdata
            clear my_data;
            [my_data(:,:) N_sweeps N_chn] = import_abf(file,-1,1/Ts);
             %% Take subset of data Cursor
%             xChan = [1 4]; %Channels to extract
%              xTime = [1.385 1.401]; %First Spike
%              xTime = [1.425 1.433]; 
%              xTime = [1.552 1.558]; 
              xTime = [.071 .080]; 
%                       xTime = [-1 -1];
             xChan = [1 3]
             if ~any(xTime == -1)
                 my_data = subsetAxonData(my_data,N_chn, xChan,xTime,Ts);
                 N_chn = size(xChan,2);
             end
             if true %VIEW DATA
                 figure(1) %Plot first and last 2 SWEEPS of  CHN 1, 2
                 subplot(1,2,1)
                 plot(my_data(:,[2 2+N_chn end+1-N_chn end+1-2*N_chn]))
                 subplot(1,2,2)
                 plot(my_data(:,[N_chn+1 (2*N_chn+1) end-N_chn end]))
             end


             N_samples = size(my_data,1);
             lastfileloaded = file;
        end
        % ---------------------------------------------------------------------
        % ----
        Intra_Chn = 1;
        Extra_Chn = 2; %index in xChan

         threshold = -30;
        index_beg_spike = FindSpikes(my_data(:,1+Intra_Chn:N_chn:end),threshold,Ts,1);
        if isempty(index_beg_spike)
            error('index_beg_spike is EMPTY.  FindSpikes found nothing');
        end
        [tempXspikes data skipped misalign_ind]= xSpikesIN(index_beg_spike,my_data,Ts,N_chn);
        %filename must have 19 columns
        filename = [blanks(19-length(filename)) filename];
        %             processedFiles = [processedFiles; filename];
        index_file = index_file +1;
        processedFiles{index_file,1} = filename; %java.lang.String(filename);
        Nspikes_in_File(i) = size(tempXspikes,2);

        extracted_spikes = [extracted_spikes tempXspikes];
        if isempty(extracted_spikes)
            error('extracted_spikes is EMPTY.  xSpikesIN extracted nothing');
        end

    end
end
processedFiles

if ~exist('EE_coord')
    error('EE_coord ENTER one then Save');
end
save(strcat(savefileheader,savefilename),'Ts','processedFiles','Nspikes_in_File','misalign_ind','extracted_spikes','EE_coord', '-mat')
savefilename

run('PlotSpikes001')
%save mean with name

%****OBJECTIVE online analysis.
