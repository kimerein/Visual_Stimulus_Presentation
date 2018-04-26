readdirheader = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Data Analysis\mdata\';
% readdirheader = 'E:\My documents\Scanziani Lab\Data Analysis\mdata\';

readfilename = []; %% for all files
%% all MAT  in directory
temp = strcat(readdirheader,'Copy*.mat');
filedata = struct2cell(dir(temp));

if isempty(filedata)
    stemp = sprintf('No Files found at: %s\nCheck path',temp)
    error(stemp);
end
sselectData = cell(0); selectData = [];selectDataNum = [];selectTs=[];
for(i = 1:size(filedata,2))  %% for all MAT files
    readfilename =   filedata{1,i} ;
    fileIN = strcat(readdirheader,readfilename);

load(fileIN,...
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