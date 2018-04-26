function [index_beg_spike] = FindSpikes2(data, threshold,Ts,b_intra)
%% Function takes in Intra/Extracellular data and find spikes that 
%% are larger(or smaller for Extra) then threshold and within the 
%% peakwidth defined by maxpeakwidth minpeakwidth
%%
%% Returns index of the beginning of the spike in data
%%
%% Differs from FindSpike in that it does not find index_end_spike, or
%% check max width of peak just finds threshold crossing.
if b_intra
    data_GT = int16((data > threshold)) ;
else
    data_GT = int16((data < threshold)) ;
end

% FIND spikes
% temp = sparse([diff(data_GT); zeros(1,size(data_GT,2))]); %insert a row because diff is small then the matrix it comes from
temp = [diff(data_GT); zeros(1,size(data_GT,2),'int16')]; %insert a row because diff is small then the matrix it comes from
% clear data_GT;

index_beg_spike = find(temp == 1);

