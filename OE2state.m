% Corey Marcus
% UT Austin, Aerospace Engineering, NEAR research group
% Script converts orbital elements to position and velocity in ECI frame

clear
close all
clc
format compact

mu = 398600; %km3/s2
a = 7226; %km
e = .02;

%angles in degrees
i = 12;
Ohm = 10;
w = 31;
theta = 10;

%calculate norm of position (norm of r)
nr = a*(1-e^2)/(1+e*cosd(theta));

%calculate angular velocity
h = sqrt(mu*nr*(1+e*cosd(theta)));

%calculate r in perifocal reference frame
r_pqw = h^2/mu/(1+e*cosd(theta))*[cosd(theta) sind(theta) 0]';

%calculate v in the perifocal reference frame
v_pqw = mu/h*[-sind(theta) e+cosd(theta) 0]';

%generate 313 rotation matrix from perifocal to geocentric vector rotation
DCM = angle2dcm(pi*Ohm/180, pi*i/180, pi*w/180, 'ZXZ')';

%translate vectors into ECI frame
r_ijk = DCM*r_pqw
v_ijk = DCM*v_pqw
