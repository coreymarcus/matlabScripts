function out = SE3Transform(se3,x)
%SE3Transform finds a transformation for se3 state consisting of
%[translation; axis-angle]

% % Extract parameters
% t = se3(1:3);
% omega = se3(4:6);
% 
% % Extract parameters from axis angle
% theta = sqrt(omega(1)^2 + omega(2)^2 + omega(3)^2 + 1E-6);
% e = omega/theta;
% 
% % find r_C2
% out = cos(theta)*(x) + sin(theta)*cross(e,x) + (1 - cos(theta))*dot(e,x)*e + t;

% Determine size of x
N = size(x,2);

% Find SE3 matrix representation
mat = SE3mat(se3);

% Transform
out = mat(1:3,1:4)*[x; ones(1,N)];

end