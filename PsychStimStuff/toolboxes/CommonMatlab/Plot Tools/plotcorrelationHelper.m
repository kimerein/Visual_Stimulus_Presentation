function time = plotcorrelationHelper(x,dt)
% function time = plotcorrelationHelper(x,dt)
% creates a time axis centered on zero to plot correlations
time = (([1:x]-((x-1)/2 +1)).*dt.*1000)';
