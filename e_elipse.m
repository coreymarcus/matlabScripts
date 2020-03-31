function [] = e_elipse(handle,P,sigma,color)
%Plot Error Elipses on figure
figure(handle)
hold on

% Calculate the eigenvectors and eigenvalues
[eigenvec, eigenval ] = eig(P);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, ~] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

%create theta
theta = 0:.01:2*pi;

%create semi-major and minor axes
a=sigma*sqrt(largest_eigenval);
b=sigma*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta );
ellipse_y_r  = b*sin( theta );

%Define a rotation matrix
R = [ cos(angle) sin(angle); -sin(angle) cos(angle) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Draw the error ellipse
plot(r_ellipse(:,1),r_ellipse(:,2),color)

%make a legend
label = strcat(string(sigma), ' \sigma Error Elipse');
legend('Data',label,'Location','best')

end

