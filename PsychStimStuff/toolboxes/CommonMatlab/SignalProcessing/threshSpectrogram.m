function [ind m_psd_frange pw_th]= threshSpectrogram(S,t,f,SDth,frange,minWindow,movingwin)
% function threshSpectrogram(S,t,f,SDth,frange,minWindow,movingwin)
% DESCRIPTION:  takes spectral data and applies a threshold (SDth) in power within
% frequency range (frange) and returns s
%
% INPUT :   
%    Spectragram data from (mtspectramc.m)
%         S - power (matrix with 2 or 3 (multiple sweeps) dimensions) X x Y
%                                                                     x Z
%         t  - time label for S (vector)
%         f  - freq label for S (vector)
%         movingwin  - [window windowstep] used to generate S
%
%    Thresholding parameters
%         SDth = st dev above or below (negative) mean to apply threshold.
%                 e.g. -0.5
%         frange = frequency range to apply threshold in (same units as f) 
%                   e.g. [20 60]
%         minWindow = minimum time window which must be above threshold to
%         be included (same units as t)
%                  e.g. 0.1
% OUTPUT :
%         ind : N by 3 matrix containing, 
%              COL 1: sweep number (i.e index of 3rd dimension of S)
%              COL 2: start time of S above threshold (NOTE: TIME not index             index)
%              COL 3: stop time of S above threshold
%
%       m_psd_frange:  mean power in frange each time t, and sweep (matrix
%       (X x Z)
%       pw_th:  value of power at SDth (scalar)
% USAGE: Specific function for selecting regions of LFP to analyzed based
% on power in frequency band.
% BA061308

ind_frange = [find(f>frange(1),1,'first') find(f<frange(2),1,'last')];
m_psd_frange = squeeze(mean(S(:,ind_frange(1):ind_frange(2),:),2))';
pw_th = mean(m_psd_frange(1:end))+ SDth*std(m_psd_frange(1:end)); % get threshold as XSD from mean
% pw_th = 1.0231e-011;
% APPLY threshold power and min time above threshold
minWindow = minWindow/(t(2)-t(1));% convert to ind of time in S
ind = thresholdplusinterval(m_psd_frange,pw_th,round(minWindow));
if ~isempty(ind)
    ind(:,1:3) = [ind(:,1) (ind(:,2:3)-1)*movingwin(2)+ movingwin(1)];% convert to time
end
