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

R = eye(3);

for ii = 1:3
    if order(ii) == 1
        R = [1 0 0;
            0 cos(angles(ii)) sin(angles(ii));
            0 -sin(angles(ii)) cos(angles(ii))]*R;
    elseif order(ii) == 2
        R = [cos(angles(ii)) 0 -sin(angles(ii));
            0 1 0;
            sin(angles(ii)) 0 cos(angles(ii))]*R;
    elseif order(ii) == 3
        R = [cos(angles(ii)) sin(angles(ii)) 0;
            -sin(angles(ii)) cos(angles(ii)) 0;
            0 0 1]*R;
    else
        disp('ERROR: Invalid Order')
        return
    end
end

end

