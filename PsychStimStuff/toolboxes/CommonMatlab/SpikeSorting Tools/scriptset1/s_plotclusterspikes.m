% plot cluster data
% % s_plotclusterspikes

% use with s_SpikeSort...m
%   after s_kmeanskluster.m
%%
%% makes n_units clusters of kmeans klusters
%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
T = cluster(Z,'maxclust',n_units) ;
toc


NN = 4; % subplots per figure
NP = 5 % plots per cluster;
MAXJ = size(unit,2);
figString = 'ClusterSpikes';
Nwindows = ceil(size(unit,2)/NN)
unit_Nspikes = zeros(1,n_units);
figIND = 300;
str = 'clu';

s_plotLusterSpikes

