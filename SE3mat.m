function [mat] = SE3mat(se3)
%SE3mat constructs the 4x4 matrix representation for se3 state consisting 
% of [translation; axis-angle]
mat_bottom = [0 0 0 0];
mat_R = CrossProductMat(se3(4:6));
mat_t = se3(1:3);
mat = [mat_R, mat_t; mat_bottom];

mat = expm(mat);

end