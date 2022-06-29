function R = FoxWrightRatio(alpha, x)
%FoxWrightRatio determines ratio of two FoxWright functions as specified in
%the supplementary material for the paper

if(alpha == 1)
    % Determine the initial term at the bottom of the ratio
    Phi = 0.5*(1 + erf(-x/2));

    % Protect against NaN when Phi is close to 1
    if(1 - Phi < 1E-10)
        Phi = 1 - 1E-10;
    end

    R = x/2 + (1/(2*sqrt(pi))) * exp(-(x^2)/4) / (1-Phi);
elseif( alpha < 1)
    % In this case we are forced to evaluate the ratio directly
    R = FoxWright((alpha + 1)/2,x,250)/FoxWright(alpha/2,x,250);
else
    % Run a recursion
    R = x/2 + (alpha-1)/2*(1/FoxWrightRatio(alpha - 1,x));
end

end