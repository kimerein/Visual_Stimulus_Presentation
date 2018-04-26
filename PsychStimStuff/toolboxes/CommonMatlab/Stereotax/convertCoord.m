function out = convertCoord(coord,theta,b_E1toE2)
% function converts between different XY coordinate systems that are offet
% by angle theta (Z is assumed to be parrelel
%
% Note direction is set for LN on Bass rig
% BA 032408
if nargin <3 || b_E1toE2==0
    out(:,1) = coord(:,2) .*sind(theta) + coord(:,1).*cosd(theta) .*-1;
    out(:,2) = coord(:,2).*cosd(theta) .*-1 + coord(:,1)*sind(theta).*-1;
    out(:,3) = coord(:,3);
else
    out(:,1) = -coord(:,2) .*sind(theta) - coord(:,1).*cosd(theta) ;
    out(:,2) = -coord(:,2).*cosd(theta) + coord(:,1)*sind(theta);
    out(:,3) = coord(:,3);
end
    
% NOTE moving in negative LN-Y axis, moves towards wall (E1)
%                positive LN_Y' axis moves towards wall (E2)
