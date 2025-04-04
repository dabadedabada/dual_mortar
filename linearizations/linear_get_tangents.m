function [delt_t1, delt_t2] = linear_get_tangents(cont_face, averaged_normals, delt__averaged_normals, t1)

nN = cont_face.info.nN;% Number of Nodes on the face
delt_t1 = zeros(nN, 3);
delt_t2 = zeros(nN, 3);

for i=1:nN
  v = [1, 0, 0];
  if norm(cross(v, averaged_normals(i,:))) < 1e-10
      v = [0, 1, 0];
  end
  delt_v = [0,0,0]; % just for checking
  t1_hat = v-dot(averaged_normals(i,:),v)*averaged_normals(i,:);
  t2_hat = cross(averaged_normals(i,:),t1(i,:));
  delt_t1_hat = delt_v - (dot(averaged_normals(i,:),delt_v)+dot(delt__averaged_normals(i,:), v))*averaged_normals(i,:)+dot(averaged_normals(i,:),v)*delt__averaged_normals(i,:);
  delt_t1(i,:) = delt_t1_hat/norm(t1_hat)-(dot(t1_hat,delt_t1_hat)*t1_hat)/norm(t1_hat)^3;
  delt_t2_hat = cross(delt__averaged_normals(i,:), t1(i,:))+cross(averaged_normals(i,:), delt_t1(i,:));
  delt_t2(i,:) = delt_t2_hat/norm(t2_hat)-(dot(t2_hat,delt_t1_hat)*t2_hat)/norm(t2_hat)^3;
end