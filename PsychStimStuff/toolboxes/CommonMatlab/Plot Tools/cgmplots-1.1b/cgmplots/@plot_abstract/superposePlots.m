function a_plot = superposePlots(plots, axis_labels, title_str, command, props)

% superposePlots - Superpose multiple plots with common command onto a single axis.
%
% Usage:
% a_plot = superposePlots(plots, axis_labels, title_str, command, props)
%
% Description:
%   The plot decoration will be taken from the last plot in the list, 
% with the exception of legend labels.
%
%   Parameters:
%	plots: Array of plot_abstract or subclass objects.
%	axis_labels: Cell array of axis label strings (optional, taken from plots).
%	title_str: Plot description string (optional, taken from plots).
%	command: Plotting command to use (optional, taken from plots)
%	props: A structure with any optional properties.
%		noLegends: If exists, no legends are created.
%		
%   Returns:
%	a_plot: A plot_abstract object.
%
% See also: plot_abstract, plot_abstract/plot, plot_abstract/plotFigure
%
% $Id: superposePlots.m 983 2008-02-20 22:15:42Z cengiz $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/09/23

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if ~ exist('props')
  props = struct;
end

command = '';
data = {};
legend = {};
if length(plots) > 1
  for one_plot = plots
    if isempty(command)
      command = one_plot.command;
    else
      if ~strcmp(command, one_plot.command)
	warning('mixed-command plot, using plot_superpose instead.');
	a_plot = plot_superpose(num2cell(plots));
	return;
      end
    end
    if isempty(data)
      data = one_plot.data;
    else
      if strcmp(command, 'bar')
        % special case for bar plots, add as columns
        % keep x-axis values from the 1st plot
        data{2} = [ data{2}, one_plot.data{2} ];
      else
        % general case (e.g., plot command) add vectors as input to command.
        data = {data{:}, one_plot.data{:}};
      end
    end
    if ~isfield(props, 'noLegends')
      legend = {legend{:}, one_plot.legend{:}};
    else
      legend = {};
    end
  end
else
  data = plots.data;
  legend = plots.legend;
end

a_plot = set(plots(1), 'data', data);
a_plot = set(a_plot, 'props', mergeStructs(props, a_plot.props));
a_plot = set(a_plot, 'legend', legend);

if exist('title_str') && ~ isempty(title_str)
  a_plot = set(a_plot, 'title', title_str);
end

if exist('command') && ~ isempty(command)
  a_plot = set(a_plot, 'command', command);
end

if exist('axis_labels') && ~ isempty(axis_labels)
  a_plot = set(a_plot, 'axis_labels', axis_labels);
end
