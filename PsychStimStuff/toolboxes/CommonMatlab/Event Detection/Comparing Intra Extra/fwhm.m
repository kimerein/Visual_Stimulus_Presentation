function [dfwhm hm xind]= fwhm(datain)
%function [dfwhm hm xind] = fwhm(datain)
% calculate full width at half max from datain  
% OUTPUT: xind - index of half max

%MAKE everything positive
offset = 0;
if any(datain<0)
    offset = min(datain);
    datain = datain-offset;
end
tmin = min(datain);
tmax = max(datain);
tmax_ind = find(datain == max(datain));

hm = (tmax - tmin)/2;
temp = abs(datain -hm);
%Find 2 indices closest to hm BEFORE peak and intrapolate
temp1 = min(temp(1:tmax_ind));
temp1_ind = find(temp(1:tmax_ind) == temp1);
temp1_ind = temp1_ind(1);
% 2nd nearest point
temp(temp1_ind) = tmax;
temp2 = min(temp(1:tmax_ind));
temp2_ind = find(temp(1:tmax_ind) == temp2);
temp2_ind = temp2_ind(1);

if(temp1_ind == min(temp2_ind,temp1_ind))
    xbeforepeak =  temp1/(temp1 + temp2);
    xbeforepeak = temp1_ind + xbeforepeak;
else
    xbeforepeak =  temp2/(temp1 + temp2);
    xbeforepeak = temp2_ind + xbeforepeak;
end

%Find 2 indices closest to hm AFTER peak and intrapolate
temp1 = min(temp(tmax_ind:end));
temp1_ind = find(temp(tmax_ind:end) == temp1);
temp1_ind = temp1_ind(1)+ tmax_ind-1;
temp(temp1_ind) = tmax ;
% 2nd nearest point
temp2 = min(temp(tmax_ind:end));
temp2_ind = find(temp(tmax_ind:end) == temp2);
temp2_ind = temp2_ind(1) + tmax_ind -1;

if(temp1_ind == min(temp2_ind,temp1_ind))
    xafterpeak =  temp1/(temp1 + temp2);
    xafterpeak = temp1_ind + xafterpeak;
else
    xafterpeak =  temp2/(temp1 + temp2);
    xafterpeak = temp2_ind + xafterpeak;
end

dfwhm =  xafterpeak - xbeforepeak;
hm = hm - offset;


if false%DEBUG
figure(998)
plot(datain);
line([xbeforepeak xafterpeak],[hm hm],'Color','r');
end

xind = [xbeforepeak xafterpeak];