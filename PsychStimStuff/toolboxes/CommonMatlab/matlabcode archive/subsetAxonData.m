function [subdata]  = subsetAxonData(axondata,N_Chn, xChn,xTime,Ts)
%% Input data from axon abf file Rows samples, Columns Sweeps (multiplexed
%% channels), Col_1 is time
%% N_Chn = Number of channels multiplexed in data
%%
%% xChn = Row Vector of channels to extract
%% xTime = [Beg End] in time of Samples to extract
%% Ts = Sample frequency

if all(xTime == -1)
    xTime(1) = 1; xTime(2) = size(axondata,1)/Ts;
end

sweeps = (size(axondata,2)-1)/N_Chn;
subdata = zeros((xTime(2)-xTime(1))*Ts+1,size(xChn,2)*sweeps+1); %allocate matrix 
subdata(:,1) = axondata(xTime(1)*Ts:xTime(2)*Ts,1); %timedata
for(ii = 1:size(xChn,2))
subdata(:,ii+1:size(xChn,2):end) = axondata(xTime(1)*Ts:xTime(2)*Ts,xChn(ii) + 1:N_Chn:end);
end