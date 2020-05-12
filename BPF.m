function [xHatOut] = BPF(sys, y, xHat, Params)
%BPF - Bootstrap Particle Filter
%   Provides an iteration of a bootstrap particle filter given a system,
%   measurement, and turning parameters. Source is Zanetti's notes on BPFs 
%   TO BE ZERO
%
% Inputs
% sys - system structure with the following elements
%   f = nonlinear state propagation function handle
%   h = function hangle mapping state to measurement
%   Q = Covariance for process noise
%   R = Covariance for measurement noise
% y = the measurement at time t = k
% xHat = the state estimate structure at time t = k-1, has the following
%   elements
%   est - [Nx1] the MMSE state estimate
%   P - [NxN] the sample MMSE estimate covariance
%   pMat - [Nxn] the matrix of particles
%   n - the number of particles
%   w - [1xN] the weight of each particle
% Params - parameter structure with the following parameters
%    Fill out if needed
% 
% Outputs
% xHatOut = the state estimate structure at time t = k

%% Setup

%local variables
n = xHat.n;
Q = sys.Q;
R = sys.R;
gaussEval = @(x, mu, P) 1/sqrt((2*pi)^length(x) * det(P))*exp(-.5*(x - mu)'*P^(-1)*(x - mu));


%% Resampling

%get all the weights
wMat = xHat.w;

%split particles
split = sysresample(wMat);

%duplicate particles
xHat.pMat = xHat.pMat(:,split);

%reassign weights
xHat.w = 1/n*ones(1,n);

%% Propagation

for ii = 1:n
    
    %get particle
    xHat_minus = xHat.pMat(:,ii);
    
    % Draw Particles From Importance Distrubution (Bootstrap - p(x_k|x_k-1)
    mu = sys.f(xHat_minus);
    xHat.pMat(:,ii) = mvnrnd(mu, Q)';
end

%% Update

for ii = 1:n
    
    % get particle and weight
    w_minus = xHat.w(ii);
    xHat_minus = xHat.pMat(:,ii);
    
    %anticipated measurement
    h = sys.h(xHat_minus);
    
    % update the weight
    xHat.w(ii) = w_minus*gaussEval(y, h, R);
    
end


% Renormalize the weights
xHat.w = xHat.w/sum(xHat.w);


%find the MMSE Estimates
xHat.est = sum(xHat.w.*xHat.pMat,2);

%Covariance
xHat.P = zeros(length(xHat.est));
for kk = 1:n
    xHat.P = xHat.P + xHat.w(ii)*(xHat.pMat(:,ii) - xHat.est)*(xHat.pMat(:,ii) - xHat.est)';
end

%% Output

xHatOut = xHat;
end

