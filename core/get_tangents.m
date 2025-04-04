function [t1,t2] = get_tangents(cont_face, averaged_normals)

nN = cont_face.info.nN; % Number of Nodes on the face
t1 = zeros(nN, 3);
t2 = zeros(nN, 3);

for i=1:nN
  v = [1, 0, 0];
  if norm(cross(v, averaged_normals(i,:))) < 1e-10
      v = [0, 1, 0];
  end
  t1(i,:) = v-dot(averaged_normals(i,:),v)*averaged_normals(i,:)/norm(v-dot(averaged_normals(i,:),v)*averaged_normals(i,:));
  t2(i,:) = cross(averaged_normals(i,:), t1(i,:))/norm(cross(averaged_normals(i,:), t1(i,:)));
end