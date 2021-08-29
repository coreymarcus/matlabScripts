clear
close all
clc

format long

% Create covariance
P = eye(3);
P(2,2) = 4;
P(3,3) = 5;
P(1,2) = 1;
P(2,3) = 1.5;
P = P + P';

% Plane z - value
zPlane = 2;
pt = [0 0 zPlane]';

% create symbolic function
syms x y z 'real'
p = gaussEval([x y z]',[0 0 0]',P);

% Integrate
symInt = int(int(int(p,x,-inf,inf),y,-inf,inf),z,zPlane, inf);
numInt = vpa(symInt,10)

% % Solve eigenvalues of P
% [V, D] = eig(P);
% 
% % Find standard deviation distance
% A = zeros(3);
% A(:,1) = V(:,1)*sqrt(D(1,1));
% A(:,2) = V(:,2)*sqrt(D(2,2));
% A(:,3) = V(:,3)*sqrt(D(3,3));
% stddevdist = norm(A\pt);
% 
% % Volume of spherical segment
% r = 1; %nominal radius
% V1 = 2*pi*r^2*(r - stddevdist)/3;
% 
% % Volume of sphere
% V2 = 4*pi*r^3/3;
% 
% % Volume of area above plane
% V3 = V1 - pi*(r^2 - stddevdist^2)*stddevdist/3;
% 
% % Fractional probability
% probPred = V3/V2

% Vector along which integration is performed
v = [0 0 1]';

% Find variance along vector
var2 = P(3,3);

% Find CDF at plane
cdfEval2 = 0.5*(1 + erf(zPlane/sqrt(2*var2)));

% Pdf approx
int2 = 1 - cdfEval2