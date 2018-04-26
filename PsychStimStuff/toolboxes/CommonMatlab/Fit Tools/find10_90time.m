function time10_90 = find10_90time(data,dt,dwin)
% function time10_90 = find10_90time(data,dt,dwin)
% time 10-90 of peak to peak amp in units of dt
% INPUT:
% data -  vector or matrix containing data to fit. If matrix then new piece of data in each row
% dt - sample interval (optional)
% dwin - window in data to use for fit (samples) (optional)
% OUTPUT:
%  time10_90= units of dt
%  BA090808
bdebug = 1;
if nargin<2 || isempty(dt)
    dt = 1;
end
if nargin<3 || isempty(dwin)
    dwin = [1 size(data,2)];
end
time10_90 = [];

 FITAMPRANGE = [.1 .9]; % fit 10%-90%
 xtime= [1:size(data,2)]*dt;
nrows = size(data,1);
for i = 1:nrows
    temp = data(i,dwin(1):dwin(2)); 
    %find 
    drange = range(temp); 
    dmax = max(temp);    dmin = min(temp);
    dmax_ind = findvalue(temp,dmax);
    dmin_ind = findvalue(temp,dmin);
       ind(1) = findvalue(temp(min(dmin_ind,dmax_ind):max(dmin_ind,dmax_ind)),dmax -  drange*(1-FITAMPRANGE(2)))+ dwin(1)+min(dmin_ind,dmax_ind)-2;
       ind(2) = findvalue(temp(min(dmin_ind,dmax_ind):max(dmin_ind,dmax_ind)),dmax - drange*(1-FITAMPRANGE(1))) + dwin(1)+min(dmin_ind,dmax_ind)-2;
       time10_90(i) = range(ind)*dt; % in units of dt
       
           if bdebug
            figure(999);clf
            plot(xtime-xtime(dwin(1)),data(i,:));hold on;
            line([max(ind)*dt max(ind)*dt]-dwin(1)*dt,[dmax dmin],'color','r');
            line([min(ind)*dt min(ind)*dt]-dwin(1)*dt,[dmax dmin],'color','g');
            title(num2str(time10_90(i)))
  
        end
end
