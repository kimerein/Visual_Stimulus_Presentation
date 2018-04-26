function width = findpeakwidth(data,dt,dwin,fracpeak)
% width (in time) of peak at fracpeak of peak-peak amp in units of dt
% INPUT:
% data -  vector or matrix containing data to fit. If matrix then new piece of data in each row
%  NOTE: data should only contain 1 peak/trough with only 2 values at fracpeak.
%  multiple other times that have values similar to fracpeakvalue will
%  cause erros
% dt - sample interval (optional)
% dwin - window in data to use for fit (samples) (optional)
% fracpeak - fraction of peak-peak amp (0<=fracpeak<=1) to find width (optional) default= 0.5
% OUTPUT:
%  width= units of dt
%  BA090808
bdebug = 1;
if nargin<2 || isempty(dt)
    dt = 1;
end
if nargin<3 || isempty(dwin)
    dwin = [1 size(data,2)];
end
if nargin<4 || isempty(fracpeak)
    fracpeak =.5;
end
width = [];

 xtime= [1:size(data,2)]*dt;
nrows = size(data,1);
for i = 1:nrows
    temp = data(i,dwin(1):dwin(2)); 

    if sum(abs(diff(diff(temp))))>3; warning('multiple trough/peaks may cause error in computed width'); end;

    %find 
    drange = range(temp); 
    dmax = max(temp);    dmin = min(temp);
    dmax_ind = findvalue(temp,dmax);
    dmin_ind = findvalue(temp,dmin);
       ind(1) = findvalue(temp(1:max(dmin_ind,dmax_ind)),dmax -  drange*fracpeak)+ dwin(1)-1;
       ind(2) = findvalue(temp(max(dmin_ind,dmax_ind):end),dmax - drange*fracpeak) + dwin(1)+max(dmin_ind,dmax_ind) -2 ;
       width(i) = range(ind)*dt; % in units of dt
       
           if bdebug
            figure;clf
            plot(xtime-xtime(dwin(1)),data(i,:));hold on;
            line([max(ind)*dt min(ind)*dt]-dwin(1)*dt,[dmax -  drange*fracpeak dmax -  drange*fracpeak],'color','r');
%             line([min(ind)*dt min(ind)*dt]-dwin(1)*dt,[dmax dmin],'color','g');
            title(['Width:' num2str(width(i))]);
        end
end