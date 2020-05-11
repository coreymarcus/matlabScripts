function [xHatOut, PhatOut] = GQKF(sys, y, xHat, Phat, Params)
% GQKF - Gaussian Quatrature Kalman Filter - provides an iteration of an
% Gauss-Hermite Quatrature Kalman Filter given a system, measurement, and
% turning parameters. Sourced from Advanced Estimation Course Notes and "I.
% Arasaratnam, S. Haykin, and R. J. Elliot. Discrete-time nonlinear
% filtering algorithms using gauss-hermite quadrature. Proceedings of the
% IEEE, 95(5):953–977, May 2007."
% 
% WARNING: FILTER IS CURRENTLY SPECIFIED FOR M = 2 AND NUMBER OF STATES =
% 4, ADDITIONAL WORK IS NEEDED TO GENERALIZE IT. SPECIFICALLY WE NEED AN
% ALGORITHM TO PERMUTATE ALL OF THE EIGENVALUES CORRECTLY
% 
% Inputs
% sys - system structure with the following elements
%   f = nonlinear state propagation function handle
%   h = nonlinear measurement function handle
%   Q = process noise covariance
%   R = measurement noise covariance
% y = the measurement at time t = k
% xHat = the state estimate at time t = k-1
% Phat = the state estimate covariance at time t = k-1
% Params - parameter structure with the following parameters
%   m = order of the gauss hermite quadrature rule, currently supported m =
%       2, 3, 4, 5
%
% Outputs
% xHatOut = the state estimate at time t = k
% PhatOut = the state estimate covariance at time t = k

%% Setup
n = length(xHat); %size of state
l = length(y); %size of measurement

%extract some locals
m = Params.m;
Q = sys.Q;
R = sys.R;

%number of quadrature points
Nq = m^n;

%% Propagation

%decompose P
S = chol(Phat,'lower');

%create J for state
v = sqrt(1:(m-1));
J = diag(v,1) + diag(v,-1);

%get eigenvalues and vectors of J
[V, D] = eig(J);

%get the weights and eigenvalues
w = abs(V(1,:)).^2;
lamda = diag(D);

%permutate to get the eigenvalue permutations
% L = 1:m;
% Y = [];
% switch m
%     case 2
%         Y = allcomb(L,L);
%     case 3
%         Y = allcomb(L,L,L);
%     case 4
%         Y = allcomb(L,L,L,L);
%     case 5
%         Y = allcomb(L,L,L,L,L);
% end
% 
% %reorder
% Y = Y';
% 
% if(isempty(Y))
%     disp('ERROR: Invalid Order Number')
% end

Y = ones(4,16);
Y(2,5:8) = 2;
Y(3,[3:4 7:8]) = 2;
Y(4,2:2:8) = 2;
Y(:,9:16) = Y(:,1:8);
Y(1,9:16) = 2;

%calculate and propagate the quadrature points, and the new estimate
quadPtsX = zeros(n,Nq);
quadWeights = zeros(1,Nq);
xHatProp = zeros(n,1);
for ii = 1:Nq
    targCol = Y(:,ii);
    quadWeights(ii) = prod(w(targCol));
    quadPtsX(:,ii) = sys.f(S*lamda(targCol) + xHat);
    xHatProp = xHatProp + quadWeights(ii)*quadPtsX(:,ii);
end

%update estimate
xHat = xHatProp;

%calculate the new estimate covariance
Phat = Q;
for ii = 1:Nq
    Phat = Phat + quadWeights(ii)*(quadPtsX(:,ii)-xHat)*(quadPtsX(:,ii)-xHat)';
end


%% Update

%decompose P
S = chol(Phat,'lower');

%calculate and propagate the quadrature points, and the new estimate
quadPtsX = zeros(n,Nq);
quadPtsY = zeros(l,Nq);
quadWeights = zeros(1,Nq);
yHat = zeros(l,1);
for ii = 1:Nq
    targCol = Y(:,ii);
    quadWeights(ii) = prod(w(targCol));
    quadPtsX(:,ii) = S*lamda(targCol) + xHat;
    quadPtsY(:,ii) = sys.h(quadPtsX(:,ii));
    yHat = yHat + quadWeights(ii)*quadPtsY(:,ii);
end

%evaluate the inovation and cross covariance matricies
Pyy = R - yHat*yHat';
Pxy = -xHat*yHat';
for ii = 1:Nq
    Pyy = Pyy + quadWeights(ii)*quadPtsY(:,ii)*quadPtsY(:,ii)';
    Pxy = Pxy + quadWeights(ii)*quadPtsX(:,ii)*quadPtsY(:,ii)';
end

%Kalman Gain
K = Pxy/Pyy;

%state update
xHat = xHat + K*(y - yHat);

%covariance update
Phat = Phat - K*Pyy*K';


%% Outout
xHatOut = xHat;
PhatOut = Phat;

end

