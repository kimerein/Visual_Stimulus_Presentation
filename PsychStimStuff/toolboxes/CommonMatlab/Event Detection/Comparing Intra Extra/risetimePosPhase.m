function [risetime] = risetimePosPhase(datain)
% function [dfwhm hm] = fwhmPosPhase(datain)
% find fwhm of positive phase of Extra
%
tmin = min(datain);
tmin_ind = find(datain == tmin);
tmax = max(datain(tmin_ind:end));
tmax_ind = find(datain(tmin_ind:end) == tmax);
tmax_ind = tmax_ind + tmin_ind -1;
temp1_ind = find(abs(datain(tmin_ind:tmax_ind)) == min(abs(datain(tmin_ind:tmax_ind)))) +tmin_ind-1
temp1 = datain(temp1_ind);
if (abs(datain(temp1_ind+1)) < abs(datain(temp1_ind-1)))
    temp2_ind = temp1_ind+1;
    temp2 = datain(temp2_ind);
else
    temp2_ind = temp1_ind-1;
    temp2 = datain(temp2_ind);
end  

if(temp1_ind == min(temp2_ind,temp1_ind))
    zerox =  temp1/(temp1 + temp2);
    zerox = temp1_ind + zerox;
else
    zerox =  temp2/(temp1 + temp2);
    zerox = temp2_ind + zerox;
end

risetime = tmax_ind - zerox;
% %find zerocrossing after negative peak and before postive peak
% temp = min(abs(datain(tmin_ind:tmax_ind)));
% temp_ind = min(find(abs(datain(tmin_ind:end))== temp(1)));
% [dfwhm hm] = fwhm(datain(temp_ind+tmin_ind-1:end));
if true %DEBUG
figure(997)
plot(datain);
line([zerox tmax_ind],[0 0]);
end