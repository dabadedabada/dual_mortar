function Phatje = linear_get_nodal_element_normals(cont_face, xi_id, ele_id)
% popp 3d (A3) or popp diss. (A.3)
% linearization of non-unit normals on nodes for each element by using two element tangent vectors
% returns matrix Phatje for single node j an element e
% delt_nhatje = Phatje*delt_d
% xi_index is index of isoparam element node, ele_index is index of element

nN_ele = cont_face.info.nN_ele; % Number of Nodes on a single element
nN = cont_face.info.nN;
fe = cont_face.fe;
Phatje = zeros(3, nN*3);
ele_nod_coo = cont_face.coo(cont_face.nod(ele_id,:),:); % coordinates of nodes of elements
deriv = fe.dNdxi(fe.xi(xi_id,:)');
for k=1:nN_ele
  % popp 3d (A3) or popp diss. (A.3)
  temp_Phatje = eye(3)*(deriv(k,1)+deriv(k,2))+skew_mat(ele_nod_coo(k,:))*(deriv(k,1)-deriv(k,2));
  Phatje(:,cont_face.nod(ele_id,k)*3-2:cont_face.nod(ele_id,k)*3) = temp_Phatje;
end