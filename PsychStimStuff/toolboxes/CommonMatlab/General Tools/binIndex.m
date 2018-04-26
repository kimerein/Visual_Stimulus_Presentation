function [ind binCenter data] = binIndex(data,binsedges)
% function [ind binCenter] = binIndex(data,binsedges)
%  INPUT: data vector that was used to determine which bin each element falls into
%     binsedges: monotonic increasing vector containing the edges of each bin
%  OUTPUT: ind vector of length N where length(data) = N
%            contains the bin which each element belongs
%          binCenter: center of each bin
%
% BA070508
ind = zeros(size(data)-1);
binCenter = zeros(size(binsedges)-1);

for i=2:length(binsedges)
    ind(find((data <= binsedges(i))&data>binsedges(i-1)))=i-1;
    binCenter(i-1) = (binsedges(i)+binsedges(i-1))/2;
end

 ind(find(data>binsedges(end)))=-1; % larger than largest bin
 ind(find(data<binsedges(1)))=-2; % smaller than smallest bin
