function [cross_corr] = spiketime_xcorr(Xdata,Ydata,nSperBin,range)
%function [cross_corr] = spiketime_xcorr(Xdata,Ydata,nSperBin,range)
% 
% Creates a raster plot with nbins for Xdata and Ydata and xcorrs them over
% the range
%
% Xdata and Ydata contain spiketimes
% 
if isequal(Xdata,Ydata)  % AUTOCORR
        % Bin Data
        binned_spikesX = bin2(Xdata(:,1),nSperBin);
        %         plot((1:size(binned_spikes,2))*(binsize),binned_spikes(1,:))  % X axis time in Se
        cross_corr = xcorrmanual(binned_spikesX,binned_spikesX,range,1);
else
        % Bin Data
        binned_spikesX = bin2(Xdata(:,1),nSperBin);
        binned_spikesY = bin2(Ydata(:,1),nSperBin);
        %         plot((1:size(binned_spikes,2))*(binsize),binned_spikes(1,:))  % X axis time in Se
        cross_corr = xcorrmanual(binned_spikesX,binned_spikesY,range,1);
end