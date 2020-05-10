%complex step derivative
function deriv = CX(x,fun,h)
%assumes fun output is scalar or vector valued
%find dimensions
n = length(x);
test = fun(x);
N = length(test);

deriv = zeros(N,n);
for ii = 1:n
    dx = zeros(n,1);
    dx(ii) = h*1i;
    funEval = fun(x + dx);
    deriv(:,ii) = imag(funEval)/h;
end

end