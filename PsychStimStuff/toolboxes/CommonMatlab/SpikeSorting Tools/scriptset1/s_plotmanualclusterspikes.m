% plot cluster data
% % s_plotmanualclusterspikes

% use with s_SpikeSort...m
%   after s_kmeanskluster.m
%%
%% makes n_units clusters of kmeans klusters

% %MANUAL select klusters that make up 1 unit
% unit = cell(n_units,1);
unit{1} = [1 6]; % A
unit{2} = [4]; % A
% unit{2} = [19 27 63 55 62]%B
% unit{3} = [15 4 2 28]%C
% unit{4} = [1]%D
% unit{5} = [3]%E
% unit{6} = [6]%F
% unit{7} = [7]%G
% unit{8} = [8]%H
% unit{9} = [9]%I
% unit{10} = [10]%J
% unit{11} = [12]%K
% unit{12} = [13]%L
% unit{13} = [14]%M
% unit{14} = [39 52 54 48 44 50 29 47 41 22 25 36 51 53 5 21]%N
% % % unit{16} = []%O
% % % unit{15} = [58 62]%P
% % %
% %
NN = 4; % subplots per figure
NP = 5 % plots per cluster;
MAXJ = size(unit,2);
figString = 'ManualClusterSpikes';
Nwindows = ceil(size(unit,2)/NN)
unit_Nspikes = zeros(1,n_units);
figIND = 400;
str = 'manclu';

s_plotLusterSpikes