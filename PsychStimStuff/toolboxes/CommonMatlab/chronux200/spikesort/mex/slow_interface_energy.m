%
% Author:  Dan Hill
%   This is a temporary function to replace the MEX function
%   "fast_interface_energy".  
%

function output = slow_interface_energy( dists,scale )

    if nargin ~= 2, error('Function takes two inputs.'); end
      
    output = sum( exp( dists(:) / -scale ) );
    
  
