function xHat = initGSF(n, x0, P0, Np)
%INITGSF Initializes a Gaussian Sum Filter based on a gaussian distribution
% Inputs
% n - number of components desired
% x0 - initial distribution mean
% P0 - initial distribution covariance
% Np - number of points to be drawn from the initial distribution
% 
% Outputs
% xHat - an initial estimate structure with the following components
%   est - [Nx1] the state estimate
%   P - [NxN] the estimate covariance
%   compMu - [NxL] the mean of each component
%   compP - [NxNxL] the covariance of each component
%   n - the number of components
%   w - [1xN] the weight of each component

%draw points from initial distribution
m = length(x0);
X = zeros(Np, m);
for ii = 1:Np
    X(ii,:) = mvnrnd(x0, P0);
end

%fit gaussian mixture distribution
GM = fitgmdist(X,n);

%build xHat
xHat.w = GM.ComponentProportion;
xHat.compMu = GM.mu';
xHat.compP = GM.Sigma;
xHat.n = n;

%calculate mean and covariance
l = length(x0);
xHat.est = zeros(l,1);
xHat.P = zeros(l);
for ii = 1:n
    xHat.est = xHat.est + xHat.w(ii)*xHat.compMu(:,ii);
    xHat.P = xHat.P + xHat.w(ii)*(xHat.compP(:,:,ii) + xHat.compMu(:,ii)*xHat.compMu(:,ii)' - xHat.est*xHat.est');
end

end

