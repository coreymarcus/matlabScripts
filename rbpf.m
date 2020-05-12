function [xHatOut, xMMSE_l, xMMSE_n, P_lOut, P_nOut] = rbpf(sys, y, xHat, Params)
%rbpf - Rao Blackwellised Particle Filter - provides an iteration of an
%   RBPF given a system, measurement, and turning parameters. Uses a
%   bootstrap importance distribution. Source is Zanetti's notes on RBPFs
%   CURRENTLY f_l ASSUMED TO BE ZERO
%
% Inputs
% sys - system structure with the following elements
%   f_n = nonlinear state propagation function handle
%   A_n = function handle to generate matrix mapping linear state to new
%       nonlinear state
%   B_n = function handle to generate matrix mapping noise to nonlinear
%       state
%   f_l = linear state propagation function handle
%   A_l = function handle to generate matrix mapping linear state to new
%       linear state
%   B_l = function handle to generate matrix mapping noise to linear
%       state
%   h = function hangle mapping nonlinear state to measurement
%   C = function handle creating matrix to map linear state to measurement
%   D = function handle creating matrix to map noise to measurement
%   Pnu_n = Covariance for nonlinear propagation noise
%   Pnu_l = Covariance for linear propagation noise
%   Peta = Covariance for measurement noise
%   N_n = dim of nonlinear state
%   N_l = dim of linear state
% y = the measurement at time t = k
% xHat = the state estimate at time t = k-1, cell array where each cell is a
%   structure corresponding to a particle with the following fields
%   w - the particle weight
%   xHat_l - an estimate of the linear states
%   P_l - the covariance of the linear states
%   xHat_n - an estimate of the nonlinear states
% Params - parameter structure with the following parameters
%   Npart = Number of particles to be used for the nonlinear states
%
% Outputs
% xHatOut = the state estimate at time t = k
% xMMSE_l - a MMSE estimate of the linear states
% xMMSE_n - a MMSE estimate of the nonlinear states
% P_lOut - the covariance of the linear states
% P_nOut - the sample covariance of the nonlinear states

%% Setup

%local variables
Npart = Params.Npart;
gaussEval = @(x, mu, P) 1/sqrt((2*pi)^length(x) * det(P))*exp(-.5*(x - mu)'*P^(-1)*(x - mu));

%check to verify f_l = 0
if(sum(abs(sys.f_l(zeros(sys.N_n,1)))) ~= 0)
    disp('ERROR: Nonzero f_l!')
    return
end

%% Resampling

%get all the weights
wMat = zeros(Npart,1);
for ii = 1:Npart
    wMat(ii) = xHat{ii}.w;
end

%split particles
split = sysresample(wMat);

%duplicate particles
xHat = xHat(split);

%reassign weights
for ii = 1:Npart
    xHat{ii}.w = 1/Npart;
end

%% Propagation

for ii = 1:Npart
    
    %local variables
    xHat_n_minus = xHat{ii}.xHat_n;
    xHat_l_minus = xHat{ii}.xHat_l;
    P_l_minus = xHat{ii}.P_l;
    A_l = sys.A_l(xHat_n_minus);
    A_n = sys.A_n(xHat_n_minus);
    B_n = sys.B_n(xHat_n_minus);
    B_l = sys.B_l(xHat_n_minus);
    
    % Draw Particles From Importance Distrubution (Bootstrap - p(x_k|x_k-1)
    mu = sys.f_n(xHat_n_minus) + A_n*xHat_l_minus;
    xHat{ii}.xHat_n = mvnrnd(mu, A_n*P_l_minus*A_n' + B_n*sys.Pnu_n*B_n')';
    
    % update the mean and covariance of all the KFs conditioned on the new
    % particles
    W = A_n*P_l_minus*A_n' + B_n*sys.Pnu_n*B_n';
    K = P_l_minus*A_n'/W;
    xHat{ii}.xHat_l = xHat_l_minus + K*(xHat{ii}.xHat_n - mu);
    xHat{ii}.P_l = P_l_minus - K*W*K';
    
    % propagate the mean and covariance of the KFs
    xHat{ii}.xHat_l = A_l*xHat{ii}.xHat_l;
    xHat{ii}.P_l = A_l*xHat{ii}.P_l*A_l' + B_l*sys.Pnu_l*B_l';
    
    
end

%% Update

%track weights for normalization
wTrack = 0;

for ii = 1:Npart
    
    % Local variables
    w_minus = xHat{ii}.w;
    xHat_n = xHat{ii}.xHat_n;
    xHat_l = xHat{ii}.xHat_l;
    R = sys.Peta;
    h_n = sys.h(xHat_n);
    C = sys.C(xHat_n);
    D = sys.D(xHat_n);
    P_l = xHat{ii}.P_l;
    
    % Evauluate the gaussian
    p = gaussEval(y, h_n + C*xHat_l, C*P_l*C' + D*R*D');
    
    % update the weights
    xHat{ii}.w = w_minus*p;
    wTrack = wTrack + xHat{ii}.w;
    
    % Update the KFs with the measurement
    W = C*P_l*C' + D*R*D';
    K = P_l*C'/W;
    xHat{ii}.xHat_l = xHat_l + K*(y - h_n - C*xHat_l);
    xHat{ii}.P_l = P_l - K*W*K';
    
end

if(wTrack == 0)
    disp('Error: Unstable Weight Normalization!')
end

% Renormalize the weights
for ii = 1:Npart
    xHat{ii}.w = xHat{ii}.w/wTrack;
end

%find the MMSE Estimates
xMMSE_l = zeros(sys.N_l,1);
xMMSE_n = zeros(sys.N_n,1);

for ii = 1:Npart
    xMMSE_l = xMMSE_l + xHat{ii}.w*xHat{ii}.xHat_l;
    xMMSE_n = xMMSE_n + xHat{ii}.w*xHat{ii}.xHat_n;
end

%Covariance
P_lOut = zeros(sys.N_l);
P_nOut = zeros(sys.N_n);
for kk = 1:Npart
    P_nOut = P_nOut + xHat{kk}.w*(xHat{kk}.xHat_n - xMMSE_n)*(xHat{kk}.xHat_n - xMMSE_n)';
    P_lOut = P_lOut + xHat{kk}.w*(xHat{kk}.P_l + xHat{kk}.xHat_l*xHat{kk}.xHat_l' - ...
        xMMSE_l*xMMSE_l');
end

%% Output

xHatOut = xHat;
end

