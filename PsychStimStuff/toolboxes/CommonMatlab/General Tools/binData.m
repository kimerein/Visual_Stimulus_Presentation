%BVA 111005
% creates a rasterdata given spike times
% creates bindata (binning data into specified windows

function [rasterdata bindata] = binData(spiketime,binwindow)
Ts = 20e-6; % assume 50kHz sampling

%% Raster is 1 when there is a spike
for iii = 1:size(spiketime,1)
    rasterdata(spiketime(iii)) = uint8(1); % keep sweeps
end
%add time column (Formatting for Binning)
% rasterdata = [[1:uint32(samplePerSweep*N_sweeps)]' rasterdata];
rasterdata = [[1:max(spiketime)]' rasterdata];

%% Bin Data
binData = binning(rasterdata,binwindow*Ts,0,size(rasterdata,1))*binwindow*Ts; %number of spikes per bin
binData = binData(:,1:2); %other columns useless
if ((sum(binData(:,2)) ~= size(spiketime,1))) %sanity check # of spikes (should be same as index_beg_spike)
    error('spike in binData ~= number of extacted spikes');
end
