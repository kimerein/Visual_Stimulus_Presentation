clear all
close all

% colororder2 =['r-';'g-';'b-';'c-';'m-';'k-'];
colororder =['r.';'g.';'b.';'c.';'m.';'k.'];
readdirheader = 'C:\Documents and Settings\Bassam\My Documents\Scanziani Lab\Data Analysis\mdata\';
% readdirheader = 'E:\My documents\Scanziani Lab\Data Analysis\mdata\';

readfilename = []; %% for all files
temp = strcat(readdirheader,'*type*.mat');
filedata = struct2cell(dir(temp));

if isempty(filedata)
    stemp = sprintf('No Files found at: %s\nCheck path',temp)
    error(stemp);
end

init = 0;
for(j = 1:size(filedata,2))
    j
    readfilename =   filedata{1,j} ;
    fileIN = strcat(readdirheader,readfilename);

    load(fileIN,...
        'processedFiles','Ts','Nspikes_in_File','extracted_spikes', 'misalign_ind', 'EE_coord','-mat');

    % mEE = zeros(size(Nspikes_in_File),size(extracted_spikes,1));
    clear temp11;
    for i = 1: size(Nspikes_in_File,2)
        ii = i + init;
        start_ind = sum(Nspikes_in_File(1,1:i-1));
        % correct for dc offset (bit of a hack)
        mEE(ii,:) = (mean(extracted_spikes(:,start_ind+1:start_ind+Nspikes_in_File(1,i),3),2)'- mean(mean(extracted_spikes(1:20,start_ind+1:start_ind+Nspikes_in_File(1,i),3),1)));
        distance(ii,1) = sqrt(EE_coord(i,2)^2+ EE_coord(i,3)^2 + EE_coord(i,4)^2);
       temp11(i) = min(mEE(ii,:));
        %     distance(i,1) = sqrt(EE_coord(indMin,2)^2+ EE_coord(indMin,3)^2 + EE_coord(indMin,4)^2);
        %     distance(i,1) = sqrt(EE_coord(indMin,2)^2+ EE_coord(indMin,3)^2 + EE_coord(indMin,4)^2);

    end


    amp(init+1:init+size(Nspikes_in_File,2)) = temp11./min(temp11);
    figure(1)
    plot(mEE');
    pause;
    figure(2);
    plot(distance(init+1:init+size(Nspikes_in_File,2)),amp(init+1:init+size(Nspikes_in_File,2)),colororder(j,:),'MarkerSize',15);
    hold on;
    ylim([0 1.1])
    
    NpositionPerFile(j) = size(Nspikes_in_File,2);
    init = init+size(Nspikes_in_File,2);

end
% title(fileIN(1+max(strfind(fileIN,'\')):end-4),'Interpreter','none')
% title('O-2_280905_typepyr_B02.mat')
xlabel('Distance um');
ylabel('Amplitude uV');

% normalize to max
