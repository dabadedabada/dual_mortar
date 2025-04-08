function Ttilde = linear_project_onto_plane(coo, nod, x0, n0, N0, C0, sl_or_m, nN_m)
% Popp 3d (A17)
% nod is as integer - index of the row in the cont_face.coo
% coo is a vector 1x3 of coordinates of point

S = zeros(size(N0));
if (sl_or_m == 1)
  index = nN_m*3 + nod*3-2;
  S(:,index:index+2) = eye(3);
end
if (sl_or_m == 2)
  index = nod*3-2;
  S(:,index:index+2) = eye(3);
end

Ttilde = S - (n0'*n0)*(S-C0) + ((n0*(coo-x0)')*eye(3)-(coo-x0)'*n0)*N0;

