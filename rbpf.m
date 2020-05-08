function [xHatOut, xMMSE_l, xMMSE_n] = rbpf(sys, y, xHat, Params)
%rbpf - Rao Blackwellised Particle Filter - provides an iteration of an
%   RBPF given a system, measurement, and turning parameters. Source is
%   Zanetti's notes on RBPFs
%   CURRENTLY A_n AND f_l ARE ASSUMEED TO BE ZERO
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
% xHat - the state estimate FILL IN MORE DETAIL
% Params - parameter structure with the following parameters
%   Npart = Number of particles to be used for the nonlinear states
%
% Outputs
% xHatOut = the state estimate at time t = k, stacked as [x_l' x_n']'

%% Setup

%local variables
Npart = Params.Npart;
gaussEval = @(x, mu, P) 1/sqrt((2*pi)^length(x) * det(P))*exp(-.5*(x - mu)'*P^(-1)*(x - mu));

%check to verify A_n and f_l are zero
if(sys.A_n ~= 0)
    disp('ERROR: Non-zero A_n')
    return
end

if(sys.f_l ~= 0)
    disp('ERROR: Non-zero f_l')
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
    B_n = sys.B_n(xHat_n_minus);
    B_l = sys.B_l(xHat_n_minus);
    
    % Draw Particles From Importance Distrubution (Bootstrap - p(x_k|x_k-1)
    mu = sys.f_n(xHat_n_minus);
    xHat{ii}.xHat_n = mvnrnd(mu, B_n*sys.Pnu_n*B_n')';
    
    % update the mean and covariance of all the KFs conditioned on the new
    % particles (Don't need to do this because A_n = 0)
    
    % propagate the mean and covariance of the KFs (this is unecessary for
    % SLAM because map is static)
    xHat{ii}.xHat_l = A_l*xHat_l_minus;
    xHat{ii}.P_l = A_l*P_l_minus*A_l' + B_l*sys.Pnu_l*B_l';
    
    
end

%% Update

%track weights for normalization
wTrack = 0;

for ii = 1:Npart
    
    % Local variables
    w_minus = xHat{ii}.w;
    xHat_n = xHat{ii}.xHat_n;
    xHat_l = xHat{ii}.xHat_l;
    %     A_l = sys.A_l(xHat_n);
    B_n = sys.B_n(xHat_n);
    %     B_l = sys.B_l(xHat_n);
    f_n = sys.f_n(xHat_n);
    Q_n = sys.Pnu_n;
    R = sys.Peta;
    h_n = sys.h(xHat_n);
    C = sys.C(xHat_n);
    D = sys.D(xHat_n);
    P_l = xHat{ii}.P_l;
    
    % Evauluate the gaussians
    %     p1 = gaussEval(xHat_n, f_n, B_n*Q_n*B_n');
    p2 = gaussEval(y, h_n + C*xHat_l, C*P_l*C' + D*R*D');
    %     p3 = gaussEval(xHat_n, f_n, B_n*Q_n*B_n');
    
    % update the weights
    %xHat{ii}.w = w_minus*p1*p2/p3;
    xHat{ii}.w = w_minus*p2; %equivalent for bootstrap
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

% % perform kalman filter updates for each particle, finding p(x_l | x_n, y)
% for ii = 1:Npart
%
%     %extract particle, build necessary matricies
%     p = xHat{ii};
%     C = sys.C(p.xHat_n);
%     D = sys.D(p.xHat_n);
%     P = p.P_l;
%     R = sys.Peta;
%     xHat_l = p.xHat_l;
%     xHat_n = p.xHat_n;
%
%     %Kalman update
%     K = P*C'/(C*P*C' + D*R*D');
%     p.P_l = (eye(sys.N_l) - K*C)*P*(eye(sys.N_l) - K*C)' + K*D*R*D'*K';
%     p.xHat_l = xHat_l + K*(y - sys.h(xHat_n) - C*xHat_l);
%
%     %reassign
%     xHat{ii} = p;
%
% end

% using the measurement, update the weights for each particle

% propagate the kalman filter for each particle

% use the propagation to sample the particles

% use the samples as a measurement to perform another kalman update

% obtain an estimate of the linear states

% obtain an estimate of the nonlinear states



xHatOut = xHat;
end

