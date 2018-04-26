function prodata = filterdata(data,dt,cutoff,type)
% function prodata = filterdata(data,dt,cutoff,type)
% Function filters data 
% INPUT:  
%     data-    matrice of data.
%     dt -    sampling (used to turn cutoff into appropriate frequency)
%     cutoff - "cut off frequency" (Hz)
%     type - 0 = low-pass, 1 = high-pass
%     
% OUTPUT:
%     prodata
%
% BA103006

if type ==0 
    stype = 'low';
elseif type ==1
    stype = 'high';
else
    error('Unknown fiter type')
end

[B,A] = butter(2,2*cutoff*dt,stype);
prodata =   filtfilt(B,A,data);