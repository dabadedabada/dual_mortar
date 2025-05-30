function R = rotation_matrix(n0)
% returns rotation matrix of Rodrigues' rotation formula

% normalize the n0 vector
n0 = n0 / norm(n0);

% Compute rotation axis (cross product with [0,0,1])
target_normal = [0; 0; 1];
rhat = cross(n0, target_normal);
r = rhat / norm(rhat);

% Compute rotation angle
theta = acos(dot(n0, target_normal));

% Check if the n0 is already aligned
if norm(rhat) < 1e-6
  R = cos(theta)*eye(3); % Identity if no rotation needed
  %R = eye(3);
else
  % Rodrigues' rotation formula
  R = r*r' + cos(theta)*(eye(3)-r*r') + sin(theta) * skew_mat(r);
end

