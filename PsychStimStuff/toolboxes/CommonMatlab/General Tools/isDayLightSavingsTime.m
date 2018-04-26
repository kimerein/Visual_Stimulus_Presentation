function [blogical startDST endDST] = isDayLightSavingsTime()
% function [blogical startDST endDST] = isDayLightSavingsTime()
% BA
DATVEC = datevec(now);
% DST starts 1st Sunday in Nov
i=1; startDST= 0;
while weekday(startDST) ~= 1 % until find sunday
    startDST = datenum([DATVEC(1) 11 i 2 0 0]);
    i= i+1;
end
% DST ends 2nd Sunday in March
i=1; j = 0;
while j<2 % until find sunday
    endDST = datenum([DATVEC(1) 3 i 2 0 0]);
    if weekday(endDST) == 1
        j = j+1;
    end
    i= i+1;
end

blogical = 0;
if now >= startDST | now < DATVEC 
    blogical = 1;
end


