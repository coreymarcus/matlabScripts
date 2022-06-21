function [se3] = SE3vec(mat)
%SE3vec converts a 4x4 se3 matrix to vector

% initialize output
se3 = zeros(6,1);

% Find elements
logmat = logm(mat);
se3(1:3) = logmat(1:3,4);
se3(4) = logmat(3,2);
se3(5) = logmat(1,3);
se3(6) = logmat(2,1);

end