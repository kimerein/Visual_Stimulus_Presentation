function nomatch = findnomatch(d,dtofind)
% function nomatch = findnomatch(data,dtofind)
%
% Find which dtofind are in d

% OUTPUT:
%         nomatch: vector with size = size(dtofind)
%         0 if no match else 1

nomatch = zeros(size(dtofind));
for i=1:length(dtofind)
    if find(d == dtofind(i)) ~= -1
        nomatch(i) = 1;
    end
end
        