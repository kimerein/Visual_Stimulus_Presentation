function [dfwhm hm] = fwhmPosPhase(datain)
% function [dfwhm hm] = fwhmPosPhase(datain)
% find fwhm of positive phase of Extra
%

tmin = min(datain);
tmin_ind = find(datain == tmin);
tmax = max(datain(tmin_ind:end));
tmax_ind = find(datain == tmax);
%find zerocrossing after negative peak and before postive peak
temp = min(abs(datain(tmin_ind:tmax_ind)));
temp_ind = min(find(abs(datain(tmin_ind:tmax_ind))== temp(1)));
[dfwhm hm] = fwhm(datain(temp_ind+tmin_ind-1:end));
