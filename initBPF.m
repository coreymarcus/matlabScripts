function [xHatOut] = initBPF(n,x0,P0)
%initBPF - initializes a Bootstrap Particle Filter with a gaussian prior
%   on the state
% 
% Inputs
% n = Number of desired particles
% x0 = mean of the initial distribution
% P0 = covariance of the initial distribution
% 
% Outputs
% xHatOut = the state estimate structure at time t = k, has the following
%   elements
%   est - [Nx1] the MMSE state estimate
%   P - [NxN] the sample MMSE estimate covariance
%   pMat - [Nxn] the matrix of particles
%   n - the number of particles
%   w - [1xN] the weight of each particle

%create the particle weights
xHatOut.w = ones(1,n)/n;
xHatOut.n = n;

%draw particles
xHatOut.pMat = zeros(length(x0),n);
xHatOut.est = zeros(length(x0),1);
for ii = 1:n
    xHatOut.pMat(:,ii) = mvnrnd(x0, P0)';
    xHatOut.est = xHatOut.est + xHatOut.pMat(:,ii)/n;
end

%Covariance
xHatOut.P = zeros(length(xHatOut.est));
for kk = 1:n
    xHatOut.P = xHatOut.P + xHatOut.w(ii)*(xHatOut.pMat(:,ii) - xHatOut.est)*(xHatOut.pMat(:,ii) - xHatOut.est)';
end


end

