function [m, b, P] = TLS(x, y, sig2_x, sig2_y)
% TLS - Total Least Squares. Performs a single dimension linear fit of a
% line y=mx+b to a set of x and y data points. Use when there is
% measurement noise in x and y, and when you are equally confident in all
% measurements. Sourced from "A tutorial on the total least squares method
% for fitting a straight line and a plane" by Munoz et al, 2014
% 
% INPUTS:
% x - an [nx1] or [1xn] vector of x measurements
% y - an [nx1] or [1xn] vector of y measurements
% sig2_x - the error variance in x measurments
% sig2_y - the error variance in y measurments
% 
% OUTPUTS:
% m - slope of best fit line
% b - y intercept of best fit line
% P - the estimate covariance matrix. I'm not currently sure how to
%     estimate this for TLS, so I've output the standard Least Squares
%     covariance matrix

%number of data points
n = length(x);

%average x and y
x_tilde = sum(x)/n;
y_tilde = sum(y)/n;

%calculate two sums
sum_num = 0;
sum_den = 0;
for ii = 1:n
    sum_num = sum_num + (x(ii) - x_tilde)*(y(ii)-y_tilde);
    sum_den = sum_den + (y(ii) - y_tilde)^2 - (x(ii) - x_tilde)^2;
end

%find phi
phi = 0.5*atan2(-2*sum_num,sum_den);

%find r
r = x_tilde*cos(phi) + y_tilde*sin(phi);

%find line parameters
m = -1/tan(phi);
b = r/sin(phi);

%estimate the covariance
H = ones(n,2);
H(:,1) = x;
P = (sig2_x+sig2_y)*eye(2)/(H'*H);

end

