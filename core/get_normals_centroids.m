function normals_storage = get_normals_centroids(cont_face)

nN = cont_face.info.nN; % Number of nodes on the face
% Get the normals by sum of nodal normals and then mapping at the center
ele_normals = get_nodal_element_normals(cont_face);
sum_normals = get_sum_normals(cont_face, ele_normals);
averaged_normals = zeros(size(sum_normals));
for i=1:nN
  averaged_normals(i,:) = sum_normals(i,:)/norm(sum_normals(i,:));
end
[t1,t2] = get_tangents(averaged_normals);

normals_storage.ele_normals = ele_normals;
normals_storage.sum_normals = sum_normals;
normals_storage.averaged_normals = averaged_normals;
normals_storage.n0 = get_n0(cont_face, averaged_normals);
normals_storage.x0 = get_centroids(cont_face);
normals_storage.t1 = t1;
normals_storage.t2 = t2;

