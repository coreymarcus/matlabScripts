function [se3out] = PropagateSE3(se3,v,omega,dt)
%PropagateSE3 propagates an se3 vector with constant velocity and rate

% Form matrix
M = SE3mat(se3);

% Extract R and t and angles
R = M(1:3,1:3);
t = -R'*M(1:3,4);
angle = rotm2eul(R','XYZ')';

% Propagate
t = t + dt*v;
angle = angle + dt*omega;
R = frameROTM(angle,[1 2 3]);

% Convert back to a matrix
Mout = [R -R*t;
    0 0 0 1];

% Convert to se3 vector
se3out = SE3vec(Mout);

end