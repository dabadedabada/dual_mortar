function n0 = get_n0(cont_face, avg_normals)
% calculates normal in centroids of elements using reference element
% parameter space (popp A.29)
% returns 3xnele non-unit n0, normals in centroids of elements ordered as mesh.bou{}.ele

nele = cont_face.info.nele;% Number of Elements on the face
fe = cont_face.fe;

shapef =  fe.N(fe.Cxi); % Shape funcs at centroid of isoparametric (reference) element
n0 = zeros(nele, 3);
for i=1:nele
  curr_el_normals = avg_normals(cont_face.nod(i,:),:);
  n0(i,:) = shapef'*curr_el_normals;
  n0(i,:) = n0(i,:)/norm(n0(i,:));
end