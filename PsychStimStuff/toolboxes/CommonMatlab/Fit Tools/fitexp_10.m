function fit_tau = fitexp_10-90(data,dt,dwin)
%function fitexp_10-90(data,dt,dwin)
% Fit data with exponential from 10-90 of peak to peak amp
% INPUT:
% data -  vector or matrix containing data to fit. If matrix then new piece of data in each row
% dt - sample interval
% dwin - window in data to use for fit (samples)
% OUTPUT:
%  fit_tau= exp fit (sec)
%  BA090808

bdebug = 1;

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
temp = data(i,:); xtime = xtime1; %end
    %find 
    drange = range(temp(dwin(1):dwin(2))); 
    dmax = max(temp(dwin(1):dwin(2)));    dmin = min(temp(dwin(1):dwin(2)));
       ind(1) = findvalue(temp(dwin(1):dwin(2)),dmax -  drange*(1-FITAMPRANGE(2))) + dwin(1);
       ind(2) = findvalue(temp(dwin(1):dwin(2)),dmax - drange*(1-FITAMPRANGE(1))) + dwin(1);
 
        [decayTau1 FittedCurve estimates sse] = fitDecaySingle(temp(min(ind):max(ind)));
        fit_tau(i) = decayTau*dt; % sec
        if bdebug
            figure(999);clf
            plot(xtime-xtime(dwin(1)),temp);hold on;
            line([max(ind)*dt max(ind)*dt]-min(ind)*dt,[dmax dmin],'color','r');
            line([min(ind)*dt min(ind)*dt]-min(ind)*dt,[dmax dmin],'color','g');
            Ha = gca();
                 fit1exp(temp(min(ind):max(ind)),dt,dt,Ha);
                 axis tight
  
        end
    
end

  
 