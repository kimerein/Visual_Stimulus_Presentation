WOI = [5 50]/1000*Ts; % window to take around event
lengtht = floor(WOI(2)-WOI(1));
evW = [];
for i = 1:size(da,1)
    b = da(i,2)-WOI(1);
    e = da(i,2)+WOI(2);
   if  (b>0)& (e< swlength)
       b = floor(b);
       e = b + lengtht;
       evW = [evW; I(da(i,1),b:e)]%event waveform
   end
end
    
