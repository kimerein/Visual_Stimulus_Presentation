function quad_ind = quartindex(data,xedge)
% function quad_ind = quartindex(data)
%  INPUT: data that was used to determine quadrent edges 
% (use quadrent function for this)
%         xedge  quadrent edges 
%  OUTPUT: quad_ind is 4 by length(data)
%           row 1 is quadrant 1.  Each element is 1 if it is part of the
%           quadrant otherwise 0

quad_ind = zeros(length(data),4,'int16');
quad_ind(:,1) = (data <= xedge(3))';
quad_ind(:,2) = (data <= xedge(2)) & (data > xedge(3));
quad_ind(:,3)  = (data <= xedge(1)) & (data > xedge(2));
quad_ind(:,4)  =  data > xedge(1);

% quad_ind = cell(4,1);
% quad_ind{1} = find( data <= xedge(3));
% quad_ind{2} = find( (data <= xedge(2)) & data > xedge(3));
% quad_ind{3} = find( (data <= xedge(1)) & data > xedge(2));
% quad_ind{4} = find( data > xedge(1));
