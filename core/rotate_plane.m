function R = rotate_plane(n0)
    % returns rotation matrix of Rodrigues' rotation formula
    % points2D = (R * points3D')';
    % pointsBack = (R' * points2D')';

    % normalize the n0 vector
    n0 = n0 / norm(n0);
    
    % Compute rotation axis (cross product with [0,0,1])
    target_normal = [0; 0; 1];
    rot_Axis = cross(n0, target_normal);
    
    % Compute rotation angle
    theta = acos(dot(n0, target_normal));

    % Check if the n0 is already aligned
    if norm(rot_Axis) < 1e-8
        R = eye(3); % Identity if no rotation needed
    else
        rot_Axis = rot_Axis / norm(rot_Axis);
        % Rodrigues' rotation formula
        K = [  0          -rot_Axis(3)  rot_Axis(2);
               rot_Axis(3)   0         -rot_Axis(1);
              -rot_Axis(2)  rot_Axis(1)  0];
        R = eye(3) + sin(theta) * K + (1 - cos(theta)) * (K * K);
    end
end
