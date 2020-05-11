function [xHatOut, PhatOut] = UKF(sys, y, xHat, Phat, Params)
%UKF - Unscented Kalman Filter - provides an iteration of an UKF given a
%system, measurement, and turning parameters.
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
%   sigPts = string with the choice in sigma point arrangment, the
%       following are supported ('2n+1')
%   k = unnormalized weight of the central point, valid for '2n+1' sig
%       points only
%
% Outputs
% xHatOut = the state estimate at time t = k
% PhatOut = the state estimate covariance at time t = k

%% Setup
n = length(xHat); %size of state
m = length(y); %size of measurement

%extract some locals
sigPts = Params.sigPts;

%build sigma points array and weight array
switch sigPts
    case '2n+1'
        Xi = zeros(n,2*n+1);
        Y = zeros(m,2*n+1);
        k = Params.k;
        w = zeros(1,2*n+1);
end

%detect size of sigma points array
L = size(Xi,2);


%Matrix square root
S = chol(Phat,'lower');

%% Propagation

%grab sigma points
switch sigPts
    case '2n+1'
        for ii = 1:n
            Xi(:,ii) = xHat - sqrt(n+k)*S(:,ii);
            Xi(:,ii + n) = xHat + sqrt(n+k)*S(:,ii);
            w(ii) = .5/(n+k);
            w(ii+n) = .5/(n+k);
        end
        
        %center point
        Xi(:,2*n+1) = xHat;
        w(2*n+1) = k/(n+k);
end

%propagate sigma points and average estimate
xHat = zeros(n,1);
for ii = 1:L
    Xi(:,ii) = sys.f(Xi(:,ii));
    xHat = xHat + w(ii)*Xi(:,ii);
end

%propagate the covariance
Phat = sys.Q;
for ii = 1:L
    Phat = Phat + w(ii)*(Xi(:,ii) - xHat)*(Xi(:,ii) - xHat)';
end

%% Update

%take the matrix square root again
S = chol(Phat,'lower');

%recalculate the sigma points based on the propagated xHat
% also evaluate Y at each sigma point
switch sigPts
    case '2n+1'
        for ii = 1:n
            Xi(:,ii) = xHat - sqrt(n+k)*S(:,ii);
            Xi(:,ii + n) = xHat + sqrt(n+k)*S(:,ii);
            Y(:,ii) = sys.h(Xi(:,ii));
            Y(:,ii + n) = sys.h(Xi(:,ii + n));
        end
        
        %center point
        Xi(:,2*n+1) = xHat;
        Y(:,2*n+1) = sys.h(Xi(:,2*n+1));
end

%get the measurement mean, and the various covariances
yHat = zeros(m,1);
for ii = 1:L
    yHat = yHat + w(ii)*Y(:,ii);
end

%get Pyy and Pxy
Pyy = sys.R;
Pxy = zeros(n,m);
for ii = 1:L
    Pyy = Pyy + w(ii)*(Y(:,ii) - yHat)*(Y(:,ii) - yHat)';
    Pxy = Pxy + w(ii)*(Xi(:,ii) - xHat)*(Y(:,ii) - yHat)';
end

%state estimate
xHatOut = xHat + (Pxy/(Pyy))*(y - yHat);

%covariance estimate
PhatOut = Phat - (Pxy/Pyy)*Pxy';

end

