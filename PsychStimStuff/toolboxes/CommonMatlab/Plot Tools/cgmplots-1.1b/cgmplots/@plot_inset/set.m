function b = set(a, attr, val)
% set - Generic method for setting object attributes.
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/10/08

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.
if ~ isa(a, 'plot_abstract')
  builtin('set', a, attr, val);
else
  try
    a.(attr) = val;
    b = a;
  catch
    errstr = lasterror;
    b = a;
    b.plot_abstract = set(a.plot_abstract, attr, val);
  end
end
