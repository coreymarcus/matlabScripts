function [fw] = FoxWright(halfalpha,x,n)
%FoxWright Calculates the special case of the Fox-Wright function as
%described in "The Modified-Half Normal Distribution: Properties..." (2021)


% Initialize
% fw = 0;
% for ii = 0:n
%     fw = fw + gamma(halfalpha + 0.5*ii)*x^ii/factorial(ii);
% end

% Try to do it faster
idxs = 0:n;
t1 = gamma(halfalpha + 0.5*idxs);
t2 = x.^idxs;
t3 = factorial(idxs);
fw = sum(t1.*(t2./t3));

end