function [T1, T2] = linear_get_tangents(cont_face, normals_storage, nod_id)
% returns the matrices for dir deriv of tangent vectors t1, t2 at node_id

ni = normals_storage.averaged_normals(nod_id,:);
vi = [1, 0, 0];
t1 = normals_storage.t1(nod_id,:);
Ni = linear_averaged_normals(cont_face,normals_storage, nod_id);
t1hat = cross(vi,ni)';
t2hat = cross(t1,ni)';
if norm(cross(vi, ni)) < 1e-10
    vi = [0, 1, 0];
end

T1hat = skew_mat(vi)*Ni;
T1 = (eye(3)/norm(t1hat)-t1hat*t1hat'*norm(t1hat)^-3)*T1hat;
T2hat = skew_mat(ni)'*T1 +skew_mat(t1)*Ni;
T2 = (eye(3)/norm(t2hat)-t2hat*t2hat'*norm(t2hat)^-3)*T2hat;


