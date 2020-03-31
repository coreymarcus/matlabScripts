% Corey Marcus
% ASE 366K - Spacecraft Dynamics
% HW #6

clear
close all
clc
format compact

%% Problem 1
disp('Problem 1')
disp(' ')

mu = 398600; %km3/s2
r = [-4743 4743 0]'; %km
v = [-5.879 -4.223 0]'; %km/s

%detirmine position and velocity 20 min later using f and g functions

%find magnitude of r and v
nr = norm(r);
nv = norm(v);

%find semi-major axis
a = -mu/(nv^2-(2*mu/nr));

%calculate magnitude of angular momentum
h = norm(cross(r,v));

%calculate eccentricity
e = sqrt(1-(h^2/mu/a));

%solve for initial eccentric anommaly
E0 = acos((1-nr/a)/e);

%solve mean motion
n = sqrt(mu/a^3);

%find out what current time since periapse is
t0 = (E0-e*sin(E0))/n;

%desired time is time since periapse plus 20 min
t = 20*60 + t0;

%define eccentric anomally function
fun = @(E) E-e*sin(E)-n*t; 

%find E at t current
E1 = fzero(fun,n*t);

%solve f and g functions
f = a/nr*((cos(E1)-e)*cos(E0)+sin(E1)*sin(E0));
g = sqrt(a^3/mu)*(sin(E1-E0)-e*(sin(E1)-sin(E0)));

%calculate new position
r1 = f*r+g*v

%calculate norm of new position
nr1 = norm(r1);

%calculate f_dot and g_dot
f_dot = -sqrt(mu*a)/nr/nr1*sin(E1-E0);
g_dot = 1 - a/nr1*(1-cos(E1-E0));

%calculate new velocity
v1 = f_dot*r+g_dot*v

%% Problem 2
% Do the same thing but with numerical integration
disp(' ')
disp('Problem 2')
disp(' ')
clear

r = [-4743 4743 0]'; %km
v = [-5.879 -4.223 0]'; %km/s

%final time
tf = 60*30;

%use diff eq solver (function at end of file
[t,y] = ode45(@grav, [0 tf], [r' v']');

r2 = y(end,1:3)'

v2 = y(end,4:6)'

%define the system of differential equations
function dxdt = grav(t,x)

mu = 398600; %km3/s2

%extract state information
r_x = x(1);
r_y = x(2);
r_z = x(3);
v_x = x(4);
v_y = x(5);
v_z = x(6);

%calculate distance from center of earth
nr = norm([r_x r_y r_z]);

%calculate accelerations
a_x = -mu*r_x/nr^3;
a_y = -mu*r_y/nr^3;
a_z = -mu*r_z/nr^3;

%generate output vector
dxdt = [v_x v_y v_z a_x a_y a_z]';

end