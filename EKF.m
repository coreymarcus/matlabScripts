function [xHatOut, PhatOut, wEval] = EKF(sys, y, xHat, Phat, Params)
%EKF - Extended Kalman Filter - provides an iteration of an EKF given a
%system, measurement, and turning parameters. Source is
% Inputs
% sys - system structure with the following elements
%   f = nonlinear state propagation function handle
%   h = nonlinear measurement function handle
%   F = jacobian of f wrt x function handle
%   H = jacobian of h wrt x function handle
%   Q = process noise covariance
%   R = measurement noise covariance
% y = the measurement at time t = k
% xHat = the state estimate at time t = k-1
% Phat = the state estimate covariance at time t = k-1
% Params - parameter structure with the following parameters
%   evalLikelihood = boolean detirmining if the likelihood of the
%       measurement is evaluated. Used when the EKF is included in the
%       Gaussian sum filter
%
% Outputs
% xHatOut = the state estimate at time t = k
% PhatOut = the state estimate covariance at time t = k

%% Setup

%extract
eval = Params.evalLikelihood;

%function to evaluate gaussians
gaussEval = @(x, mu, P) 1/sqrt((2*pi)^length(x) * det(P))*exp(-.5*(x - mu)'*P^(-1)*(x - mu));


%% Propagation

%jacobian of system dynamics
F = sys.F(xHat);

%state and covariance propagation
xBar = sys.f(xHat);
Pbar = F*Phat*F' + sys.Q;

%% Update

%evaluate measurement jacobian at xBar
H = sys.H(xBar);

%anticipated measurement
yBar = sys.h(xBar);

%innovation covariance
S = H*Pbar*H' + sys.R;

%kalman gain
K = (Pbar*H')/S;

%state estimate
xHatOut = xBar + K*(y - yBar);

%covariance estimate (Josephs Form)
I = eye(length(xHat));
PhatOut = (I - K*H)*Pbar*(I - K*H)' + K * sys.R *K';

%evaluate weight (if needed)
if(eval)
    wEval = gaussEval(y, yBar, S);
else
    wEval = [];
end

end

