function smdata = smgaussian(data,s)
% function smdata = smgaussian(data,s)
% INPUT :   data is a N x 2 array
%               column 1 is x value
%               column 2 y values
%           s is standard deviation of smoothing gaussian
% OUTPUT:   smdata is a N x 2 array
%               column 1 is x value
%               column 2 y values
%           convolved with gaussian
%% NOTE: double check that x axis is alligned corretly
bdebug = 0;

xstep = data(2,1) - data(1,1);

x = [-1*4*s:xstep:4*s]; %% x values for gaussian
g = 1/(sqrt(2.*pi).*s).*exp(-x.^2./(2*s^2));
smdata(:,2) = conv(g,data(:,2));
% align so that smdata starts at same x value as data does.
% smdata(:,1) = [((length(g)-1)*-1)*xstep:xstep:(length(smdata)-length(g))*xstep];
smdata(:,1) = [((length(g)-1)/2*-1)*xstep:xstep:(length(smdata)-length(g)/2)*xstep];

% DEBUG
if bdebug
    figure(3)
    plot(smdata(:,1),smdata(:,2))
    hold on
    plot(x,g)
end
%%%%%%%%%%%%%%%
