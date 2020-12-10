function p = gaussEval(x, mu, P) 
%gaussEval evaulates a gaussian distribution with mean mu, covariance P,
%at x

p = 1/sqrt((2*pi)^length(x) * det(P))*exp(-.5*(x - mu)'/P*(x - mu));

if(isnan(p))
    p = 0;
end

end

