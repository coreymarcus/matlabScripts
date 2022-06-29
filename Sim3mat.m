function [mat] = Sim3mat(sim3)
%Sim3mat constructs the 4x4 matrix representation for sim3 state consisting 
% of [scale; translation; rotation]
scale = sim3(1);
se3 = sim3(2:7);

% Find SE3 portion
mat = SE3mat(se3);

% Modify with scale
mat(1:3,1:3) = scale*mat(1:3,1:3);

end