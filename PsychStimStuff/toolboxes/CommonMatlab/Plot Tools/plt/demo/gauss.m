% gauss.m ---------------------------------------------------------
% This example script plots the results of combining uniform
% random variables.
% - Shows the advantage of setting the line data after the plt(..) call.
% - Note the use of the 'FigName' and 'TraceID' arguments.
% - Note the appearance of the greek letter in the x-axis label.
% - Uses the 'COLORdef' argument to select Matlab's default plotting
%   colors (typically set to use a white background for the plotting area)
% - Shows how to use the 'Options' argument to enable the x-axis cursor
%   slider (which appears just below the peak and valley finder buttons).
% - Uses the 'DIStrace' argument to initialize with some traces disabled.

% ----- Author: ----- Paul Mennen
% ----- Email:  ----- paul@mennen.org

  mxN = 10;                    % sum up to 10 uniform distributions
  dis = [0 0 0 ones(1,mxN-3)]; % initially just show the first 3 traces
  h = plt(1,ones(1,mxN),'LabelX','Standard deviation (\sigma)','LabelY','',...
           'TraceID',['Gauss'; reshape(sprintf('Sum%2d',2:mxN),5,mxN-1)'],...
           'COLORdef','default','FigName','Sum of uniform distribtions',...
           'DIStrace',dis,'xlim',[-4 4],'ylim',[-.05 1.05],'Options','S-X');
  sz = 100;          % size of each uniform distribution
  u = ones(1,sz);    % uniform distribution
  y = u;             % y will be composite distribution
  for n = 2:length(h)
    y = conv(y,u);   % convolve with next uniform distribution
    m = length(y);  mean = (m+1)/2;   sigma = sz * sqrt(n/12);
    x = ((1:m) - mean) / sigma; % change units to sigma (zero mean)
    set(h(n),'x',x,'y',y/max(y));
  end;
  set(h(1),'x',x,'y',exp(-(x.^2)/2)); % gaussian distribution

