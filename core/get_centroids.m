function centroids = get_centroids(cont_face)
% given the mesh and face (as a number) computes centroids of its
% elements 2D elements

nele = cont_face.info.nele;% Number of Elements on the face
fe = cont_face.fe;

% we use center, but the derivatives are constant (if shape funcs are linear)
shapef =  fe.N(fe.Cxi); % Shape funcs at centroid of isoparametric (reference) element
centroids = zeros(nele, 3);

for i=1:nele
  ele_nod_coo = cont_face.coo(cont_face.nod(i,:), :); % coordinates of nodes of elements
  centroids(i,:) = shapef'*ele_nod_coo;
end

