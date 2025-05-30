function Ralg = linear_rotation_matrix(n0, N0)

% Compute rotation axis (cross product with [0,0,1])
target_normal = [0; 0; 1];
rhat = cross(n0, target_normal);
r = rhat / norm(rhat);

% Compute rotation angle
Ralg = zeros(3,3,size(N0,2));
% Check if rotation was needed
if norm(rhat) < 1e-6
  Ralg(1,1,:) = N0(3,:);
  Ralg(2,2,:) = N0(3,:); 
  Ralg(3,3,:) = N0(3,:); 
  return
end

% Define the 3D matrix E (3x1x3)
E = zeros(3, 1, 3);  % Initialize a 3x1x3 matrix filled with zeros

% Define unit vectors e_x, e_y, e_z
e_x = [1; 0; 0];  % Unit vector in x-direction
e_y = [0; 1; 0];  % Unit vector in y-direction
e_z = [0; 0; 1];  % Unit vector in z-direction

% Fill the matrix E with these unit vectors in the third dimension
E(1, :, :) = e_x';  % First slice: e_x
E(2, :, :) = e_y';  % Second slice: e_y
E(3, :, :) = e_z';  % Third slice: e_z

E_skew = zeros(3,3,3);
E_skew(1,2,:) = -e_z';
E_skew(1,3,:) = e_y';
E_skew(2,1,:) = e_z';
E_skew(2,3,:) = -e_x';
E_skew(3,1,:) = -e_y';
E_skew(3,2,:) = e_x';


Rhat = (eye(3)/norm(rhat)-rhat*rhat'*norm(rhat)^(-3))*skew_mat(e_z)'*N0;

Ralg1_temp = (1-n0(3))*(pagemtimes(E,r') + pagemtimes(r, permute(E, [2, 1, 3])));
Ralg1 = pagemtimes(permute(Ralg1_temp,[1,3,2]),Rhat);
Ralg1 = permute(Ralg1,[1,3,2]);

temp = eye(3) - r*r' - skew_mat(r)*n0(3)/sqrt(1-n0(3)*n0(3));
Ralg2 = zeros(3,3,size(N0,2));
for i=1:3
  for j=1:3
    Ralg2(i,j,:) = temp(i,j)*N0(3,:);
  end
end
Ralg3_temp = pagemtimes(permute(E_skew,[1,3,2])*sqrt(1-n0(3)*n0(3)),Rhat);
Ralg3 = permute(Ralg3_temp,[1,3,2]);
Ralg = Ralg1 + Ralg2 +Ralg3;






