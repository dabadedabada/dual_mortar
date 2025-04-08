function C0 = linear_get_centroid(cont_face, ele_id)
% returns 3x3*nN matrix C0, delt_x0 = C0*delt_d
nN_ele = cont_face.info.nN_ele;
nN = cont_face.info.nN;
fe = cont_face.fe;
C0 = zeros(3,nN*3);
shapef =  fe.N(fe.Cxi); % Shape funcs at center of isoparametric (reference) element
glob_id = cont_face.nod(ele_id,:)*3-2; % vector of nN_ele node indexes of element
for k=1:nN_ele
  C0(:,glob_id(k):glob_id(k)+2) = eye(3)*shapef(k);
end



