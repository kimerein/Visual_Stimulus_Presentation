function  [fn abf_ind] = getabfindex(filepath)
%function [fn abf_ind] = getabfindex(filepath)
% extract abf filename and index of filename from path

tmp = findstr(filepath,'\');
if isempty(tmp)
    tmp =0;
end
fn = filepath(tmp(end)+1:end);
abf_ind = fn(end-7:end-4);