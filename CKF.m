function [xHatOut, PhatOut] = CKF(sys, y, xHat, Phat, Params)
% CKF - Cubature Kalman Filter - provides an iteration of an Cubature
% Kalman Filter given a system, measurement, and turning parameters.
% Sourced from Advanced Estimation Course Notes and "Ienkaran Arasaratnam
% and Simon Haykin. Cubature kalman filters. IEEE Transactions on Automatic
% Control, 54(6):1254 – 1269, June 2009. doiI: 10.1109/TAC.2009.2019800."
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
% Params - parameter to be used if necessary
%
% Outputs
% xHatOut = the state estimate at time t = k
% PhatOut = the state estimate covariance at time t = k

%% Setup
n = length(xHat); %size of state
l = length(y); %size of measurement

%extract some locals
Q = sys.Q;
R = sys.R;

%number of cubature points
Nc = 2*n;

%% Propagation

%decompose P
S = chol(Phat,'lower');

%find cubature points
xi = sqrt(Nc/2);
Xi = [xi*eye(n), -xi*eye(n)];
for ii = 1:n
    Xi(ii,ii) = xi;
    Xi(ii,ii+n) = -xi;
end

% Propagate the cubature points
cubePtsX = zeros(n,Nc);
for ii = 1:Nc
    cubePtsX(:,ii) = sys.f(xHat + S*Xi(:,ii));
end

%update estimate
xHat = sum(cubePtsX,2)./Nc;

%calculate the new estimate covariance
Phat = Q - xHat*xHat';
for ii = 1:Nc
    Phat = Phat + cubePtsX(:,ii)*cubePtsX(:,ii)'./Nc;
end


%% Update

%decompose P
S = chol(Phat,'lower');

%calculate and propagate the cubature points, and the new estimate
cubePtsX = zeros(n,Nc);
cubePtsY = zeros(l,Nc);
yHat = zeros(l,1);
for ii = 1:Nc
    cubePtsX(:,ii) = S*Xi(:,ii) + xHat;
    cubePtsY(:,ii) = sys.h(cubePtsX(:,ii));
    yHat = yHat + cubePtsY(:,ii)./Nc;
end

%evaluate the inovation and cross covariance matricies
Pyy = R - yHat*yHat';
Pxy = - xHat*yHat';
for ii = 1:Nc
    Pyy = Pyy + cubePtsY(:,ii)*cubePtsY(:,ii)'./Nc;
    Pxy = Pxy + cubePtsX(:,ii)*cubePtsY(:,ii)'./Nc;
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

