function ele_normals = get_nodal_element_normals(cont_face)
% calculates non-unit normals on nodes for each element by using two element tangent vectors
% returns nele cells of nN x 3 matrices of normals

nele = cont_face.info.nele;% Number of Elements on the face
nN_ele = cont_face.info.nN_ele;% Number of Nodes on a single element
fe = cont_face.fe;
ele_normals = cell(nele,1);

for i=1:nele
  normals_list = zeros(nN_ele,  3);
  ele_nod_coo = cont_face.coo(cont_face.nod(i,:),:); % coordinates of nodes of elements
  for j=1:nN_ele
    deriv = fe.dNdxi(fe.xi(j,:)');
    x_xi = deriv(:,1)'*ele_nod_coo;
    x_eta = deriv(:,2)'*ele_nod_coo;
    normals_list(j,:) = cross(x_eta, x_xi); % perpendicular outward vector
  end
  ele_normals{i} = normals_list;
end