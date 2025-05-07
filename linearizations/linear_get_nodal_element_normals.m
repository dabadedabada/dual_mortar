function Nhatje = linear_get_nodal_element_normals(cont_face, xi_id, ele_id)
% popp 3d (A3) or popp diss. (A.3)
% linearization of non-unit normals on nodes for each element
% returns matrix Phatje for single node j an element e
% delt_nhatje = Phatje*delt_d
% xi_in is index of isoparam element node, ele_index is index of element

nN_ele = cont_face.info.nN_ele; % Number of Nodes on a single element
nN = cont_face.info.nN;
fe = cont_face.fe;
Nhatje_xi = zeros(3, nN*3); Nhatje_eta = zeros(3, nN*3);
x_xi = zeros(3,1); x_eta = zeros(3,1);

glob_id = cont_face.nod(ele_id,:)*3-2; % vector of nN_ele node indexes of element
deriv = fe.dNdxi(fe.xi(xi_id,:)');
%x_sl_cont = reshape(cont_face.coo', nN*3, []); % reshape the slave contact nodal coordinates into vector

for k=1:nN_ele
  Nhatje_xi(:,glob_id(k):glob_id(k)+2) = eye(3)*deriv(k,1);
  Nhatje_eta(:,glob_id(k):glob_id(k)+2) = eye(3)*deriv(k,2);
  x_xi = x_xi + deriv(k,1)*cont_face.coo(cont_face.nod(ele_id,k),:)';
  x_eta = x_eta + deriv(k,2)*cont_face.coo(cont_face.nod(ele_id,k),:)';
end

Nhatje = skew_mat(x_xi)'*Nhatje_eta + skew_mat(x_eta)*Nhatje_xi;
%Nhatje = skew_mat(Nhatje_xi*x_sl_cont)'*Nhatje_eta + skew_mat(Nhatje_eta*x_sl_cont)*Nhatje_xi;
