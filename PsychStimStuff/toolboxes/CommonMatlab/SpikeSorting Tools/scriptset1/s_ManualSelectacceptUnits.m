% s_acceptUnits.m
% use in s_SpikeSort8...
%     after kmeankluster and plotcluster (or plotmanualcluster)


%Manually enter cluster to accept or join and accept
%%% ACCEPT/JOIN
clear accept_units;
accept_units{1} = {1}; %% this number is the cluster or kluster (depending on whic s_plotLusterSpike was last run)
accept_units{2} = {2};
accept_units{3} = {3};
accept_units{4} = {5};
accept_units{5} = {6};
accept_units{6} = {7};  %% removeKluster is not correct if not kluster
accept_units{7} = {8};
accept_units{8} = {9};
accept_units{9} = {10};
% accept_units{10} = {16};
% accept_units{11} = {1};
% accept_units{12} = {2};
% accept_units{13} = {6};
% accept_units{14} = {14};
% accept_units{15} = {15};
% accept_units{16} = {20};
% accept_units{17} = {21};
% accept_units{18} = {25};
% accept_units{19} = {26};
% % accept_units{18} = {15};
% % accept_units{19} = {24};
% % accept_units{20} = {18 22};
% % accept_units{21} = {25};
% % accept_units{22} = {2};
% % accept_units{22} = {2};
good_type = [1:length(GoodK)];

% accept_units{10} = {19};
%% REJECT
% reject_unit = [16:20 7 8];
reject_k=[BadK];
reject_unit = [];

sunitfilename = 'Units';