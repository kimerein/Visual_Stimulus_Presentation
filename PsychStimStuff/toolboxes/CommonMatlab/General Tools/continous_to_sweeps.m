function df = continous_to_sweeps(df,sweep_length)
% function df = continous_to_sweeps(df,sweep_length)
% converts vector into matrix of size X by sweep_length, where X =
%  ceil(length(df)/sweep_length)
%
% USAGE: for converting long contious data sets into sweeps to make
% processing less memory intensive
%
% BA
if ~isvector(df); error('df must be vector'); end
    temp = zeros(sweep_length,ceil(length(df)/sweep_length));
    temp(1:length(df)) = df;
    df = temp';
