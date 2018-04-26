function a_plot = plot_bars(mid_vals, lo_vals, hi_vals, n_vals, x_labels, y_labels, ...
			    title, axis_limits, props)

% plot_bars - Bar plot with error lines in individual axes for each variable.
%
% Usage:
% a_plot = plot_bars(mid_vals, lo_vals, hi_vals, n_vals, x_labels, y_labels, ...
%		     title, axis_limits, props)
%
% Description:
%   Subclass of plot_stack. The plot_abstract/plot command can be used to
% plot this data. Rows of *_vals will create grouped bars, columns will
% create new axes.
%
%   Parameters:
%	mid_vals: Middle points of error bars.
%	lo_vals: Low points of error bars.
%	hi_vals: High points of error bars.
%	n_vals: Number of samples used for the statistic (Optional).
%	x_labels, y_labels: Axis labels for each bar group. Must match with data columns.
%	title: Plot description.
%	axis_limits: If given, all plots contained will have these axis limits.
%	props: A structure with any optional properties.
%	  dispBarsLines: Choose between using 'bars' or 'lines' to connect the errorbars.
%	  dispErrorbars: If 1, display errorbars for lo_vals and hi_vals deviation from mid_vals 
%		     (default=1).
%	  dispNvals: If 1, display n_vals on top of each bar.
%	  groupValues: Array of within-group numeric labels, instead of just a sequence of numbers.
%	  truncateDecDigits: Truncate labels to this many decimal digits.
%	  barAxisProps: props passed to plot_abstract objects with bar commands
%		
%   Returns a structure object with the following fields:
%	plot_abstract
%
% General operations on plot_bars objects:
%   plot_bars	- Construct a new plot_bars object.
%
% Additional methods:
%	See methods('plot_bars')
%
% See also: plot_abstract, plot_abstract/plot
%
% $Id: plot_bars.m 988 2008-03-07 21:34:39Z cengiz $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/10/07

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if nargin == 0 % Called with no params
   a_plot = struct;
   a_plot = class(a_plot, 'plot_bars', plot_stack);
 elseif isa(mid_vals, 'plot_bars') % copy constructor?
   a_plot = mid_vals;
 else
   if ~ exist('props')
     props.rotateXLabel = 45; % Degrees
     %props.XTickLabel = 1;
   end
   
   group_locs = 1:size(mid_vals, 1);
   if isfield(props, 'groupValues') && ...
	 ~(isfield(props, 'XTickLabel') && isempty(props.XTickLabel))
     if isfield(props, 'truncateDecDigits')
       dig_exp = 10^props.truncateDecDigits;
       props.groupValues = round(dig_exp * props.groupValues) / dig_exp;
     end
     props.XTickLabel = num2cell(props.groupValues);
   end

   if ~ exist('axis_limits')
     axis_limits = []; % Degrees
   end

   a_plot = struct;

   num_plots = size(mid_vals, 2);
   plots = cell(1, num_plots);
   if isfield(props, 'barAxisProps')
       bar_axis_props = mergeStructs(props.barAxisProps, props);
   else
       bar_axis_props = props;
   end
   % Loop for each item and create a horizontal stack of plots
   for plot_num=1:num_plots

     if ~isfield(props, 'dispBarsLines') || strcmp(props.dispBarsLines, 'bars')
       plot_components = ...
           {plot_abstract({group_locs, mid_vals(:,plot_num)}, ...
                          {x_labels{plot_num}, y_labels{plot_num}}, '', ...
                          {}, 'bar', bar_axis_props)};
       linestyle = 'none';
     elseif strcmp(props.dispBarsLines, 'lines')
       % Enforce errorbar display then
       props.dispErrorbars = 1;
       plot_components = {};
       linestyle = '-';
     else
       error([ 'Optional argument dispBarsLines has unknown value: ' ...
               props.dispBarsLines ]);
     end

     if ~isfield(props, 'dispErrorbars') || props.dispErrorbars == 1
       plot_components = ...
	   {plot_components{:}, ...
	    plot_abstract({group_locs, mid_vals(:,plot_num), ...
			   lo_vals(:,plot_num), hi_vals(:,plot_num), 'LineStyle', linestyle}, ... % '+'
			  {x_labels{plot_num}, y_labels{plot_num}}, '', {}, 'errorbar', props)};
     end

     if ~isfield(props, 'dispNvals') || props.dispNvals == 1
       plot_components = ...
	   {plot_components{:}, ...
	    plot_abstract({group_locs, ...
			   mid_vals(:,plot_num) + hi_vals(:,plot_num), ...
                           cellfun(@(x) sprintf('n=%d', x), num2cell(n_vals(:,plot_num)), ...
                                   'UniformOutput', false), ...
			   'HorizontalAlignment', 'center', ...
			   'VerticalAlignment', 'bottom', ...
                           'FontSize', 8}, ...
			  {x_labels{plot_num}, y_labels{plot_num}}, ...
			  '', {}, 'text', props)};
     end

     plots{plot_num} = plot_superpose(plot_components, {}, '', ...
				      mergeStructs(props, struct('noLegends', 1)));
   end

   a_plot = class(a_plot, 'plot_bars', ...
		  plot_stack(plots, axis_limits, 'x', title, props));
end

% cellstr(strcat('n=', num2str(n_vals(:,plot_num))))