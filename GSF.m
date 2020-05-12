function [xHatOut] = GSF(sys, y, xHat, Params)
%GSF - Gaussian Sum Filter - provides an iteration of an GSF given a
%system, measurement, and turning parameters. Currently, each gaussian
%component is updated with an EKF
% Inputs
% sys - system structure with the following elements
%   f = nonlinear state propagation function handle
%   h = nonlinear measurement function handle
%   Q = process noise covariance
%   R = measurement noise covariance
% y = the measurement at time t = k
% xHat - an initial estimate structure with the following components
%   est - [Nx1] the state estimate
%   P - [NxN] the estimate covariance
%   compMu - [NxL] the mean of each component
%   compP - [NxNxL] the covariance of each component
%   n - the number of components
%   w - [1xN] the weight of each component
% Params - parameter structure with the following parameters
%
% Outputs
% xHatOut - the state estimate structure at time t = k


%% Local Variables
%number of components
n = xHat.n;

%EKF parameters
ParamsEKF.evalLikelihood = true;

%length of x
l = length(xHat.est);

%% Estimate
% simply run an EKF iteration for every component of xHat
pw = zeros(1,n); %vector of previous weights times liklihood
for ii = 1:n
    
    %extract
    targMu = xHat.compMu(:,ii);
    targP = xHat.compP(:,:,ii);
    targW = xHat.w(ii);
    
    %filter
    [xHatOut, PhatOut, wEval] = EKF(sys, y, targMu, targP, ParamsEKF);
    
    %update
    xHat.compMu(:,ii) = xHatOut;
    xHat.compP(:,:,ii) = PhatOut;
    
    %weight
    pw(ii) = wEval*targW;
    
end



%% Weight Update

%normalize
pw = pw/sum(pw);
for ii = 1:n
   
    %update weight
    xHat.w(ii) = pw(ii);  
    
end

%% Calculate total estimate
xHat.est = zeros(l,1);
xHat.P = zeros(l);
for ii = 1:n
    xHat.est = xHat.est + xHat.w(ii)*xHat.compMu(:,ii);
    xHat.P = xHat.P + xHat.w(ii)*(xHat.compP(:,:,ii) + xHat.compMu(:,ii)*xHat.compMu(:,ii)' - xHat.est*xHat.est');
end


%% Output
xHatOut = xHat;

end