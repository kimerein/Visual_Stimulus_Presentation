function fid = plotev(data,evtime_ind,dt,ylvl,fid)
% plotev(data,evtime_ind,dt,thresh,fid)
%
% Plots data and detected events with y-value = ylvl
%
% BA030207
if nargin<3
    dt =1;
end
xtime = [1:length(data)]*dt;

if nargin<4
    thresh = mean(data) + 3*std(data);
end

if nargin<5
    figure;
else
    figure(fid);
    clf
end

plot(xtime,data)
hold on
plot(evtime_ind.*dt,ones(size(evtime_ind)).*ylvl,'.r')
xlim([0 1])
