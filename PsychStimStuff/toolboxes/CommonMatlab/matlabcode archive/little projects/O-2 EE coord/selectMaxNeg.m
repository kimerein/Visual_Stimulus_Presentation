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
