% Corey Marcus
% UT Austin, Aerospace Engineering, NEAR research group
% Script converts position and velocity in an ECI frame to orbital elements

clear
close all
clc
format compact

mu = 398600.5; %km3/s2

r = [-2436.45, -2436.45, 6891.037]'; %km
v = [5.088611, -5.088611, 0]'; %km/s

%normalized r and v
nr = norm(r);
nv = norm(v);

%angular momentum vector
h = cross(r,v);

%semi-major axis [km]
a = -mu/(nv^2-2*mu/nr)

%eccentricity vector
e = (cross(v,h)-mu*r/nr)/mu;

%ecentricity
ne = norm(e)

%inclination
i = acosd(h(3)/norm(h))

%line of nodes
n = cross([0 0 1]',h);

%calculate Ohmega
if n(2) >= 0
    Ohm = acosd(n(1)/norm(n))
elseif n(2) < 0
    Ohm = 360-acosd(n(1)/norm(n))
else
    disp('ERROR!!!')
end

%normalized eccentricity and line of nodes
nhat = n/norm(n);
ehat = e/norm(e);

%calculate w
if e(3) >= 0
    w = acosd(nhat'*ehat) %Q1 or Q4
elseif e(3) < 0
    w = 360-acosd(nhat'*ehat) %Q2 or Q3
else
    disp('ERROR!!!')
end

%normalized radius
rhat = r/nr;

%calculate true anomally
if v'*r >= 0
    theta = acosd(ehat'*rhat) %Q1 or Q4
elseif v'*r < 0
    theta = 360 - acosd(ehat'*rhat) %Q2 or Q3
else
    disp('ERROR!!!')
end