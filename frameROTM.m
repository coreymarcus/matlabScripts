function [R] = frameROTM(angles,order)
%frameROTM creates a rotation matrix for frame rotations
%   angles is a vector containing the angles of rotation in radians
%
%   order is a vector which instructs on the order of rotation (ex: [3 2 1]
%   or [1 2 3])

%check for correct definition of order
if length(order) ~= 3
    disp('ERROR: Incorrect order definition')
    return
end

R_int = zeros(3,3,3);

for ii = 1:3
    if order(ii) == 1
        R_int(:,:,ii) = [1 0 0;
                        0 cos(angles(ii)) sin(angles(ii));
                        0 -sin(angles(ii)) cos(angles(ii))];
    elseif order(ii) == 2
        R_int(:,:,ii) = [cos(angles(ii)) 0 -sin(angles(ii));
                         0 1 0;
                         sin(angles(ii)) 0 cos(angles(ii))];
    elseif order(ii) == 3
        R_int(:,:,ii) = [cos(angles(ii)) sin(angles(ii)) 0;
                        -sin(angles(ii)) cos(angles(ii)) 0;
                         0 0 1];
    else
        disp('ERROR: Invalid Order')
        return
    end
end

R = R_int(:,:,3)*R_int(:,:,2)*R_int(:,:,1);

end

