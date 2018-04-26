function [spike_ind xclud_ind] = xcludeSpikes(datain,Ts)
% function spike_ind = xcludeSpikes(datain)
% gives indices of intracellular spikes that are most similar 
% datain is ONLY the intracellular data
%
% Algorithm: 
%       1)threshold is found from the max and min of the mean data
%       2)beg_spike holds when datain crosses threshold
%       3) the mode of beg_spike is extracted from historgram of beg_spike
%       4) spike that occur with WOI of mode are good

nbins = size(datain,2)/10; % bin size for histogram (to find mode)
WOI = 50e-6;  % sec window time around mode that spikes can occur
WOI_ind = WOI*Ts;


temp = mean(datain,2);
threshold  = max(temp) - (abs(min(temp))+ max(temp))/2;
temp = datain>threshold;
for(i = 1:size(datain,2))
    beg_spike(i) = min(find(temp(:,i)==1));
end

%Creat histogram get MODE
[temp xout]= hist(beg_spike,nbins); %Find most common threshold
mode = round(xout(find(temp==max(temp))));
%Find beg_spikes in WOI
temp = abs((beg_spike - mode * ones(1,size(beg_spike,2))));
spike_ind = find(temp < WOI_ind==1);
xclud_ind = find(temp < WOI_ind==0);

% DEBUGGING
% close 100 ;
% figure(100)
% plot(datain(:,spike_ind),'r')
% hold on;
% plot(datain(:,xclud_ind),'b')

% plot(datain)