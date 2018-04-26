function out = convertXYZ_E2(T_E2,theta,e2_XZangle,Ref_E1LNCoord,e1_XZangle)
% function out = convertXYZ_E2(T_E2,theta,e2_XZangle,Ref_LNE1,e1_XZangle)
% convert from manipulator E1XYZ coordinates to E2)
% optionally corrects for Reference provided in E1LNCoordinates
b_E1toE2 =1;
if  nargin >= 5
    if ~isempty(Ref_E1LNCoord)
        Ref = repmat(convertLNtoXYZ(Ref_E1LNCoord,e1_XZangle).*[1 1 -1],size(T_E2,1),1); %XYZ_E1% Calculate E1_Target (E1LN_coord)
    else
        Ref = repmat([0 0 0],size(T_E2,1),1);
    end
end
temp = convertCoord((T_E2 - Ref).*repmat([1 1 -1],size(T_E2,1),1),theta,b_E1toE2);
out = convertXYZtoLN(temp,e2_XZangle,[],-1);% convert to LN_E2
