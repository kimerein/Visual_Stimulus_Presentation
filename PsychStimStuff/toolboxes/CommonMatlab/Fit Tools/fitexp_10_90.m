function fit_tau = fitexp_10_90(data,dt,dwin)
%function fitexp_10-90(data,dt,dwin)
% Fit data with exponential from 10-90 of peak to peak amp
% INPUT:
% data -  vector or matrix containing data to fit. If matrix then new piece of data in each row
% dt - sample interval
% dwin - window in data to use for fit (samples)
% OUTPUT:
%  fit_tau= exp fit units of dt
%  BA090808
%%% DO: exp2fit Works better than fitDecaySingle should update

bdebug = 1;
fit_tau = [];

 FITAMPRANGE = [.1 .9]; % fit 10%-90%
%  resample = dt/10;  % for bspline
     xtime1 = [1:size(data,2)]*dt;
nrows = size(data,1);
for i = 1:nrows
%     if bspline % UNTESTED
%         xtime = [min(xtime1):resample:max(xtime1)];
%         temp = SPLINE(xtime1,data(i,:),[min(xtime1):resample:max(xtime1)]) 
%         
%     else;  
   temp = data(i,dwin(1):dwin(2));  xtime = xtime1; %end
    %find 
   drange = range(temp); 
    dmax = max(temp);    dmin = min(temp);
    dmax_ind = findvalue(temp,dmax);
    dmin_ind = findvalue(temp,dmin);
       ind(1) = findvalue(temp(min(dmin_ind,dmax_ind):max(dmin_ind,dmax_ind)),dmax -  drange*(1-FITAMPRANGE(2)))+ dwin(1)+min(dmin_ind,dmax_ind);
       ind(2) = findvalue(temp(min(dmin_ind,dmax_ind):max(dmin_ind,dmax_ind)),dmax - drange*(1-FITAMPRANGE(1))) + dwin(1)+min(dmin_ind,dmax_ind);
 
        [decayTau FittedCurve estimates sse] = fitDecaySingle(data(i,min(ind):max(ind)));
        fit_tau(i) = decayTau*dt; % sec
        if bdebug
            figure(999);clf
            plot(xtime-xtime(dwin(1)),data(i,:));hold on;
            line([max(ind)*dt max(ind)*dt]-min(ind)*dt,[dmax dmin],'color','r');
            line([min(ind)*dt min(ind)*dt]-min(ind)*dt,[dmax dmin],'color','g');
            Ha = gca();
                 fit1exp(temp(min(ind):max(ind)),dt,dt,Ha);
                 axis tight
  
        end
    
end

  
 