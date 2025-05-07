function [t1,t2] = get_tangents(averaged_normals)

nN = size(averaged_normals,1); % Number of Nodes on the face
t1 = zeros(nN, 3);
t2 = zeros(nN, 3);

for i=1:nN
  vi = [1, 0, 0];
  ni = averaged_normals(i,:);
  if norm(cross(vi, ni)) < 1e-10
      vi = [0, 1, 0];
  end
  t1hat = cross(vi,ni);
  t1(i,:) = t1hat/norm(t1hat);
  t2hat = cross(t1(i,:),ni);
  t2(i,:) = t2hat/norm(t2hat);
end