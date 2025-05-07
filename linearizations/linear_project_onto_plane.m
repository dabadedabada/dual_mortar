function P = linear_project_onto_plane(coo, nod_id, x0, n0, N0, C0, sl_or_m, nN_m)
% coo is 3x1 coo of node
% nod is as integer - index of the row in the cont_face.coo
% x0 is 3x1 centroid coo
% n0 is 3x1 normal
% N0 is 3x3(nNm + nNs) representing dir der of n0
% C0 is 3x3(nNm + nNs) representing dir der of x0
% sl_or_m is 1 or 2 - 1 if slave node, 2 if master node
% nN_m is number of nodes on master surface
% coo is a vector 3x1 of coordinates of point

S = zeros(size(N0));
if (sl_or_m == 1)
  con_id = nN_m*3 + nod_id*3-2;
elseif (sl_or_m == 2)
  con_id = nod_id*3-2;
end
S(:,con_id:con_id+2) = eye(3);
P = S - (n0*n0')*(S-C0) - ((n0'*(coo-x0))*eye(3)+n0*(coo-x0)')*N0;

