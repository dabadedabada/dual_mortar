function Phatje = linear_get_nodal_element_normals(cont_face, xi_id, ele_id)
% popp 3d (A3) or popp diss. (A.3)
% linearization of non-unit normals on nodes for each element
% returns matrix Phatje for single node j an element e
% delt_nhatje = Phatje*delt_d
% xi_in is index of isoparam element node, ele_index is index of element

nN_ele = cont_face.info.nN_ele; % Number of Nodes on a single element
nN = cont_face.info.nN;
fe = cont_face.fe;
Phatje_xi = zeros(3, nN*3);
Phatje_eta = zeros(3, nN*3);
glob_id = cont_face.nod(ele_id,:)*3-2; % vector of nN_ele node indexes of element
deriv = fe.dNdxi(fe.xi(xi_id,:)');
x_sl_cont = reshape(cont_face.coo', nN*3, []); % reshape the slave contact nodal coordinates into vector

for k=1:nN_ele
  Phatje_xi(:,glob_id(k):glob_id(k)+2) = eye(3)*deriv(k,1);
  Phatje_eta(:,glob_id(k):glob_id(k)+2) = eye(3)*deriv(k,2);
end

Phatje = skew_mat(Phatje_xi*x_sl_cont)'*Phatje_eta + skew_mat(Phatje_eta*x_sl_cont)*Phatje_xi;
