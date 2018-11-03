function [L] = Lmatrix(theta)
% Transformation matrix for a rotation of theta radians for the relation
% sigma_new = L*sigma
%
% L - Transformation matrix
% theta - rotation in the plane [radians]
L = [cos(theta)^2   sin(theta)^2    sin(2*theta);
     sin(theta)^2   cos(theta)^2    -sin(2*theta);
     -cos(theta)*sin(theta) cos(theta)*sin(theta)   cos(2*theta)];
end

