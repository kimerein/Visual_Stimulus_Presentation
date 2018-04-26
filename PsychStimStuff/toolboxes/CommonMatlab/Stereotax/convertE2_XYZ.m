function out = convertE2_XYZ(E2_coord,theta,e2_XZangle, Ref_E1LNCoord,e1_XZangle)
%function out = convertE2_XYZ(E2_coord,theta,e2_XZangle, Ref_E1LNCoord,e1_XZangle)
% convert from manipulator 2(E2) to E1XYZ coordinates)
% optionally corrects for Ref provided in E1LNCoordinates
out = convertCoord(convertLNtoXYZ(E2_coord,e2_XZangle,[],1),theta);
if  nargin >= 5
    if ~isempty(Ref_E1LNCoord)
        out= out+ repmat(convertLNtoXYZ(Ref_E1LNCoord,e1_XZangle),size(out,1),1);
    end
end

out = out.*repmat([1 1 -1],size(out,1),1);