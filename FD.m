%Finite Difference derivative
function deriv = FD(x,fun,h)
%assumes fun output is scalar or vector valued
%find dimensions
n = length(x);
test = fun(x);
N = length(test);

deriv = zeros(N,n);
parfor ii = 1:n
    dx = zeros(n,1);
    dx(ii) = h;
    funEvalF = fun(x + dx);
    funEvalB = fun(x - dx);
    deriv(:,ii) = (funEvalF - funEvalB)/(2*h);
end

end