function out = convertXYZ_E1(T_E1,e1_XZangle,Ref_E1LNCoord)
% function out = convertXYZ_E1(T_E2,e1_XZangle,Ref_E1LNCoord)
% convert from manipulator E1XYZ coordinates to E1)
% optionally corrects for Reference provided in E1LNCoordinates
Ref = repmat([0 0 0],size(T_E1,1),1);
if  nargin >= 3
    if ~isempty(Ref_E1LNCoord)
        Ref = repmat(convertLNtoXYZ(Ref_E1LNCoord,e1_XZangle).*[1 1 -1],size(T_E1,1),1); %XYZ_E1% Calculate E1_Target (E1LN_coord)
  end
end
temp = (T_E1 - Ref).*repmat([1 1 -1],size(T_E1,1),1);
out = convertXYZtoLN(temp,e1_XZangle,[],1);% convert to LN_E1
