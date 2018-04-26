function out = cLNcoord(xyz,z_xangle)
% function out = cLNcoord(xyz,z_xangle)
% account for the fact that z and x values read from LN are not orthogonal 
x = xyz(1)*sind(z_xangle);
y = xyz(2);
z = xyz(3) + xyz(1)*cosd(z_xangle);
out = [x y z];