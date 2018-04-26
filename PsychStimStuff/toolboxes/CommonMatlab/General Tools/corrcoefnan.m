function [r,p,rlo,rup] = corrcoefnan(x,y,varargin)
% function [r,p,rlo,rup] = corrcoefnan(x,y,varargin)
% performs pairwise deletion where data is missing i.e. there are NaNs
% WARNING: that this can introduce errors
% The result is only valid if the occurence of NaN's is uncorrelated. In
% order to avoid this pitfall, the correlation of NaN's should be checked 
% or case-wise deletion should be applied. 
%   Case-Wise deletion can be implemented 
%    ix = ~any(isnan([X,Y]))
%    [...] = CORRCOEF(X(ix,:),Y(ix,:),...) 
%
%  Correlation (non-random distribution) of NaN's can be checked with 
%       [nan_R,nan_sig]=corrcoef(X,isnan(X))
%   or  [nan_R,nan_sig]=corrcoef([X,Y],isnan([X,Y]))
%   or  [R,p,ci1,ci2] = CORRCOEF(...);
% also see 
% [15] http://www.statsoft.com/textbook/stbasic.html#Correlationsk
% others
ind = unique(find(~isnan(x) & ~isnan(y)));
x=x(ind);y=y(ind);

 [r,p,rlo,rup] =corrcoef(x,y,varargin{:});