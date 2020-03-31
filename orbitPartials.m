%Takes partials of potential energy for orbit propogation

clear
close all
clc

%define symbolics
syms mu R J2 x y z 'real'

%define U
r2 = x^2 + y^2 + z^2;
U = mu/sqrt(r2) * (1 - R^2/r2*J2*(3*z^2/(2*r2) - 1/2));

%take partials
dUdx = diff(U,x);
dUdy = diff(U,y);
dUdz = diff(U,z);

%define rddot
rddot = matlabFunction([dUdx, dUdy, dUdz]');
save('rddot','rddot');

