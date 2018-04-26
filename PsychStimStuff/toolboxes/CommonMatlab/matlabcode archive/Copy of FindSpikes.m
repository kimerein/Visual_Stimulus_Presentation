function [index_beg_spike] = FindSpikes(data, threshold,Ts,b_intra)
%% Function takes in Intra/Extracellular data and find spikes that 
%% are larger(or smaller for Extra) then threshold and within the 
%% peakwidth defined by maxpeakwidth minpeakwidth
%%
%% Returns index of the beginning of the spike in data

if b_intra
    maxpeakwidth = 2e-3; %sec
%     minpeakwidth = .05e-3;
     minpeakwidth = 0e-3;
    data_GT = sparse((data > threshold)) ;
else
%     maxpeakwidth = 1e-3; %sec 
    maxpeakwidth = 2e-3; %sec 
    minpeakwidth = 0e-3;
    data_GT = sparse((data < threshold)) ;
end

pack;
% FIND spikes
temp = sparse([diff(data_GT); zeros(1,size(data_GT,2))]); %insert a row because diff is small then the matrix it comes from
% clear data_GT;
% CHECK every spike has a beg before an end
for (i =1: size(temp,2))
    minIndex_beg = min(find(temp(:,i) == 1));
    minIndex_end = min(find(temp(:,i) == -1));
    if isempty(minIndex_beg) 
        if ~isempty(minIndex_end)
            temp(minIndex_end,i) = 0;
        end
    elseif isempty(minIndex_end) 
        if ~isempty(minIndex_beg)
            temp(minIndex_beg,i) = 0;
        end
    else
        if ((minIndex_end < minIndex_beg))  % end > beg
            temp(minIndex_end,i) = 0;
        end
        maxIndex_beg = max(find(temp(:,i) == 1));
        maxIndex_end = max(find(temp(:,i) == -1));
        if ((maxIndex_beg > maxIndex_end ))  % last beg has an end
            temp(maxIndex_beg,i) = 0;
        end
        %DEBUG
        if true
            tbeg = find(temp(:,i) == 1);
            tend = find(temp(:,i) == -1);
            if (size(tbeg,1) ~= size(tend,1))
                i
                error('beg and end not same')
            end
        end
    end
end
index_beg_spike = find(temp == 1);
index_end_spike = find(temp == -1);
clear temp;

%% CHOOSE SPIKES with min and maxpeakwidth
%% check that end > begin for each sweep
%% check that begin has an end
if (size(index_end_spike,1) < size(index_beg_spike,1) ) %CASE Rising edge with no falling edge
    index_beg_spike = index_beg_spike(1:end-1);
    error('ind_beg and ind_end not equal i.e.size(index_end_spike,1) < size(index_beg_spike,1')'
end

peakwidth = index_end_spike - index_beg_spike;
mask = (peakwidth > minpeakwidth*Ts) .* (peakwidth < maxpeakwidth*Ts); %Indices that meet criteria have value 1
index_beg_spike = nonzeros(mask .* index_beg_spike);
% clear mask;clear peakwidth;clear index_end_spike;
