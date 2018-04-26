function [width] = widthnegPhase(datain)
%function [ddecaytime] = widthnegPhase(datain)
% find width of 0xing 1 to 0xing 2  of Extra
%
%MAKE everything positive
offset = 0;
if any(datain<0)
    offset = min(datain);
    datain = datain-offset;
end
%%%FIND base of trough (10%)
tmin = min(datain);
tmin_ind = find(datain == tmin);
tmax = max(datain(1:tmin_ind)); % pretrough MAX
tmax_ind = find(datain(1:tmin_ind) == tmax);

%DEFINE Base of Trough as % 10per below pretoughMax
ten1 = tmax - (tmax - tmin)*.4
temp = abs(datain -ten1);

%Find 2 indices closest to  BEFORE trough and intrapolate
temp1 = min(temp(tmax_ind:tmin_ind));
temp1_ind = find(temp(tmax_ind:tmin_ind) == temp1);
temp1_ind = temp1_ind(1)+tmax_ind-1;
% 2nd nearest point
if (datain(temp1_ind+1) < datain(temp1_ind-1))
    temp2_ind = temp1_ind+1
    temp2 = temp(temp2_ind);
else
    temp2_ind = temp1_ind-1
    temp2 = temp(temp2_ind);
end    

if(temp1_ind == min(temp2_ind,temp1_ind))
    xbeforetrough =  temp1/(temp1 + temp2);
    xbeforetrough = temp1_ind + xbeforetrough;
else
    xbeforetrough =  temp2/(temp1 + temp2);
    xbeforetrough = temp2_ind + xbeforetrough;
end

%Find 2 indices closest to ten1 AFTER trough BEFORE peak and intrapolate

tmax2 = max(datain(tmin_ind:end));
tmax2_ind = find(datain(tmin_ind:end) == tmax2) +tmin_ind-1;

temp1 = min(temp(tmin_ind:tmax2_ind));
temp1_ind = find(temp(tmin_ind:tmax2_ind) == temp1);
temp1_ind = temp1_ind(1)+ tmin_ind-1;
temp(temp1_ind) = tmax ;
% 2nd nearest point
if (datain(temp1_ind+1) < datain(temp1_ind-1))
    temp2_ind = temp1_ind+1
    temp2 = temp(temp2_ind);
else
    temp2_ind = temp1_ind-1
    temp2 = temp(temp2_ind);
end  

if(temp1_ind == min(temp2_ind,temp1_ind))
    xaftertrough =  temp1/(temp1 + temp2);
    xaftertrough = temp1_ind + xaftertrough;
else
    xaftertrough =  temp2/(temp1 + temp2);
    xaftertrough = temp2_ind + xaftertrough;
end
width = (xaftertrough - xbeforetrough);
% Find base of trough
figure(999)
plot(datain);
line([xbeforetrough xaftertrough],[ten1 ten1]);
