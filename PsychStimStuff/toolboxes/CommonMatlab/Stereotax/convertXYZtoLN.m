function out = convertXYZtoLN(coord,XZangle,Ref,zdirn)
% converts to FROM orthogonal XYZ coordinates to LN coordinates using XZangle
% optionally then changes to coordinates relative to LN coordinates (if Ref
% is supplied)
% BA 032008
if nargin < 4
zdirn =1;
end
out(:,1) = coord(:,1)./sind(XZangle) ; 
out(:,2) = coord(:,2);
out(:,3) =  coord(:,3) + zdirn.*out(:,1).*cosd(XZangle); 

if nargin >2
    if ~isempty(Ref)
        out = out - repmat(Ref,size(coord,1),1);
    end
end
