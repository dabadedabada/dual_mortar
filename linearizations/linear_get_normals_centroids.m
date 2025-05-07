function delt_normals_storage = linear_get_normals_centroids(cont_face, delt_d, normals_storage)

delt_ele_normals = linear_get_nodal_element_normals(cont_face, delt_d); % Popp 3d (A3)
delt_sum_normals = linear_get_sum_normals(cont_face, normals_storage.ele_normals, delt_ele_normals); % Popp 3d (A2)
delt_averaged_normals = linear_averaged_normals(normals_storage.sum_normals,delt_sum_normals); % Popp 3d (A1)

averaged_normals = zeros(size(normals_storage.sum_normals));
for i=1:cont_face.info.nN
  averaged_normals(i,:) = normals_storage.sum_normals(i,:)/norm(normals_storage.sum_normals(i,:));
end

delt_n0 = linear_get_n0(cont_face, averaged_normals, delt_averaged_normals); % Popp 3d (A15)
delt_x0 = linear_get_centroid(cont_face, delt_x);
[delt_t1, delt_t2] = linear_get_tangents(cont_face, averaged_normals,...
  delt_averaged_normals, normals_storage.t1);

delt_normals_storage.averaged_normals = delt_averaged_normals;
delt_normals_storage.n0 = delt_n0;
delt_normals_storage.x0 = delt_x0;
delt_normals_storage.t1 = delt_t1;
delt_normals_storage.t2 = delt_t2;


