function xl = xlimhist(x,f)
% function xl = xlimhist(x,f)
%
% takes output of hist and returns the min and max x values where that have
% frequency >0
xl(1) = max(x(1),x(find(f>0,1,'first'))) ;
tmp = find(f>0,1,'last');
if tmp~= length(x)
    tmp = tmp+1;
end
xl(2) = x(tmp);