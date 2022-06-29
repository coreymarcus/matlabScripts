function out = Sim3Transform(sim3,x)
%SE3Transform finds a transformation for sim3 state consisting of
%[scale; translation; rotation]

% Number of points
N = size(x,2);

% Find SE3 matrix representation
mat = Sim3mat(sim3);

% Transform
x2 = mat*[x; ones(1,N)];

% Output
out = x2(1:3,:);

end