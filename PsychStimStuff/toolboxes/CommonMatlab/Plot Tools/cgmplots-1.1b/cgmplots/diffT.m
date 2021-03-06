function deriv = diffT(x, dy)

%  diffT - Estimate of first derivative using Taylor expansion.
%
% Usage:
% deriv = diffT(x, dy)
%
% Parameters:
%	x: A vector.
%	dy: The resolution of the discrete points in the vector.
%
% Returns:
% 	deriv: Estimate of the first derivative.
%
% Description:
%   dx     x(k-2) - 8 * x(k-1) + 8 * x(k+1) - x(k+2)
%  ---- = ------------------------------------------
%   dy			12 * dy
%
%   Taken from Sekerli, Del Negro, Lee and Butera. IEEE Trans. Biomed. Eng.,
%	51(9): 1665-71, 2004.
% Note: First and last two values of the deriv vector will contain boundary 
%	artifacts.
%
% $Id: diffT.m,v 1.3 2005/05/08 00:13:40 cengiz Exp $
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/11/15

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if size(x, 1) > size(x, 2)
  transposed = true(1);
  x = x';
else
  transposed = false(1);
end

x8 = 8 * x;

deriv = ...
    ([0, 0, 0, 0, x] - [0, 0, 0, x8, 0] + [0, x8, 0, 0, 0] - [x, 0, 0, 0, 0]) ./ ...
    ( 12 * dy );

% Strip off the boundaries
deriv = deriv(3:(end-2));

if transposed
  deriv = deriv';
end