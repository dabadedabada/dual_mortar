function Pj = linear_averaged_normals(cont_face, normals_storage, nod_id)
% return directional derivative matrix of averaged nodal normal
% delt_nj = Pj*delt_d_s
% Popp 3D (A1)

curr_normal = normals_storage.sum_normals(nod_id,:)';
Pj = (1/norm(curr_normal)-curr_normal*curr_normal'/norm(curr_normal)^3)...
  *linear_get_sum_normals(cont_face, nod_id, normals_storage.ele_normals);