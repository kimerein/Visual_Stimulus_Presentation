function [out time] = xcovm(x,y,range,dt)
% functions runs xcov on each column of x and y:
% i.e. x(:,1) with y(:,2)
%      x(:,N) with y(:,N)
%      
%      INPUT:
%         x, y: vector or matrix of data.
%             size(x,2)MUST equal size(y,2
%         range:   to run xcov over (samples)
%         dt: (optional) : sample time in (s)to create time in real time units
%              default = 1;
%       OUTPUT:
%           out :  where each column is the xcov of the respective x and y
%           columns
%           time : vector of lagtimes, length(time) = size(out,1)
if nargin <4
    dt =1;
end

if size(x,2)~=size(y,2)
    error('Two matrices x and y must have the same number of columns')
end
%predeclare 
out = zeros(range*2+1,size(x,2));

for i = 1:size(x,2) % do cov for each column
out(:,i) =     xcov(x(:,i),y(:,i),range,'coeff');
end
% create vector of lags
time = (([1:size(out,1)]-((size(out,1)-1)/2 +1)).*dt.*1000)';
