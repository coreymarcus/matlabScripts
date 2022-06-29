function [sim3] = Sim3vec(mat)
%Sim3vec converts a 4x4 se3 matrix to vector

% initialize output
sim3 = zeros(7,1);

% Find scale
sim3(1) = norm(mat(1:3,1));

% Convert mat to an se3 matrix
se3mat = mat;
se3mat(1:3,1:3) = se3mat(1:3,1:3)/sim3(1);

% Find se3 elements
sim3(2:7) = SE3vec(se3mat);

end