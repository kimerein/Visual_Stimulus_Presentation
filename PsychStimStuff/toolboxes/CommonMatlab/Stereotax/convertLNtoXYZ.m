function out = convertLNtoXYZ(coord,XZangle,Ref,zdirn)
% function corrects for Ref coordinates (if Ref is supplied)
% and then converts to orthogonal XYZ coordinates using XZangle
% BA 032008
if nargin >2
    if ~isempty(Ref)
    coord = coord + repmat(Ref,size(coord,1),1);
    end
end

if nargin < 4
zdirn =-1;
end
out(:,1) = coord(:,1) *sind(XZangle); 
out(:,2) = coord(:,2);
out(:,3) = coord(:,3) + zdirn.*coord(:,1) *cosd(XZangle); 
% note moving in negative LN-X axis, moves Down
%                positive LN-Z axis, moves Down

