%BVA 111005
% creates a rasterdata given spike times
% creates bindata (binning data into specified windows

function [rasterdata bindata] = binData2(spiketime,bin,Ts)
% DOES NOT WORK
%function [rasterdata bindata] = binData(spiketime,bin,Ts)
% Input: 
%       spiketime is a vector were each element is the # of samples, from
% the beginning of data acquistion, when a spike occured
%       binwindow is size of bins msec
%       Ts is samples per second
if nargin <3
    Ts = 50e-6 % assume 20kHz sampling
end
%% Raster is 1 when there is a spike
    temp = zeros(1,max(max(spiketime)),'int8');
    temp(1,unit_spiketimes) = 1;
%add time column (Formatting for Binning)
% rasterdata = [[1:uint32(samplePerSweep*N_sweeps)]' rasterdata];
rasterdata = [[1:max(spiketime)]' rasterdata];

%% Bin Data
binData = binning(rasterdata,binwindow*Ts,0,size(rasterdata,1))*binwindow*Ts; %number of spikes per bin
binData = binData(:,1:2); %other columns useless
if ((sum(binData(:,2)) ~= size(spiketime,1))) %sanity check # of spikes (should be same as index_beg_spike)
    error('spike in binData ~= number of extacted spikes');
end
