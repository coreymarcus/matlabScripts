function [xhist, thist] = RK4(t0, dt, N, x0, dxdt)

% Initialize output
Nstate = length(x0);
xhist = zeros(N,Nstate);
xhist(1,:) = x0';
thist = zeros(N,1);
thist(1) = t0;

for ii = 2:N
    
    % Find coefficients
    k1 = dxdt(thist(ii-1),xhist(ii-1,:)');
    k2 = dxdt(thist(ii-1) + dt/2,xhist(ii-1,:)' + k1*dt/2);
    k3 = dxdt(thist(ii-1) + dt/2,xhist(ii-1,:)' + k2*dt/2);
    k4 = dxdt(thist(ii-1) + dt,xhist(ii-1,:)' + dt*k3);

    % Update
    thist(ii) = thist(ii-1) + dt;
    xhist(ii,:) = (xhist(ii-1,:)' + (k1 + 2*k2 + 2*k3 + k4)*dt/6);
end

end