function Phatj = linear_get_sum_normals(cont_face, nod_id, ele_normals)
% returns matrix 3 x 3*nN Phatj, delt_nhatj = Phatj*delt_d_s
% linearization of sum element normals on given node
% Popp 3D (A2)

nN = cont_face.info.nN;         % Number of Nodes on the face
Phatj = zeros(3,3*nN);

[ele_id, xi_id] = find(cont_face.nod==nod_id);

for e = 1:length(ele_id)
  curr_normal = ele_normals{ele_id(e)}(xi_id,:)';
  Phatj = Phatj + (1/norm(curr_normal)-curr_normal*curr_normal'/norm(curr_normal)^3)...
    *linear_get_nodal_element_normals(cont_face,xi_id(e),ele_id(e));
end