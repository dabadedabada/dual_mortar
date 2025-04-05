function N0 = linear_get_n0(cont_face, normals_storage, ele_id)
% Popp 3d (A15) similar to (A1)
% return matrix N0, delt_n0 = N0*delt_d_s

nN = cont_face.info.nN;
nN_ele = cont_face.info.nN_ele;
fe = cont_face.fe;
shapef =  fe.N(fe.Cxi); % Shape funcs at center of isoparametric (reference) element
Nhat = zeros(3,nN*3);
curr_ele_nod = cont_face.nod(ele_id,:);
curr_ele_normals = normals_storage.averaged_normals(curr_ele_nod,:);
n0hat = curr_ele_normals'*shapef;
for k=1:nN_ele
  Nhat = Nhat + shapef(k)*linear_averaged_normals(cont_face,normals_storage,curr_ele_nod(k));
end
N0 = (1/norm(n0hat)-n0hat*n0hat'/norm(n0hat)^3)*Nhat;
  

