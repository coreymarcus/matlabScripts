function expmC = ExpmSkewSym(C)
%ExpmSkewSym produces the matrix exponential of a skew symetric matrix C

% Extract elements
a1 = C(3,2);
a2 = C(1,3);
a3 = C(2,1);
x = sqrt(a1^2 + a2^2 + a3^2);

% Matrix
expmC = eye(3) + sin(x)/x*C + (1 - cos(x))/x^2*C*C;
end