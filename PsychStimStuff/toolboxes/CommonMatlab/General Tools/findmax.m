function [index dataset] =  findmax(data1,data2,index1,index2)
% function [index dataset] =  findmax(data1,data2,index1,index2)
% performs max function on 2 values from 2 datasets.
% the values to compare are specified by the indices in index1 and index2.
%
% INPUT: 
%         data1, data2: matrix (N,M) but will be treated as linear array
%         index1, index2: vector holding indices in data1 and data2 to
%         compare.
%                   size(index1) must equal size(index2)
% OUTPUT:
%          index: holds the index of the larger value in either data1 or
%          data2
%          dataset: 1 if index is in data1 else 2.
%               e.g. if data1(index1(2))>=data2(index2(2))
%                   index(2) = index1(2)
%                   and dataset = 1;

if size(index1)~= size(index2)
    error('index1 and index2 not the same size')
end
index = zeros(size(index1));
dataset = index;
for i=1:length(index1)
    if data1(index1(i))>=data2(index2(i))
        index(i) = index1(i);
        dataset(i) = 1;
    else
        index(i) = index2(i);
        dataset(i)= 2;
    end
end
