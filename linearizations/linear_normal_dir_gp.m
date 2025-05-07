function Nhat = linear_normal_dir_gp(G_s, s_proj_gp, cf_sl, normals_storage, nod_id, nN_m)

fe = cf_sl.fe;
nN_ele_s = length(nod_id);
Nhat = zeros(size(G_s));
sum_normals = normals_storage.sum_normals(nod_id,:)';
dN_in_gp = fe.dNdxi(s_proj_gp);
N_in_gp = fe.N(s_proj_gp);
for k=1:nN_ele_s
  avg_normal = sum_normals(:,k)/norm(sum_normals(:,k));
  Nk = [zeros(3,3*nN_m),linear_averaged_normals(cf_sl, normals_storage, nod_id(k))];
  Nhat = Nhat + avg_normal*dN_in_gp(k,:)*G_s(1:2,:)+N_in_gp(k)*Nk;
end
