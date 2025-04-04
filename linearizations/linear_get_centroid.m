function delt_centroid = linear_get_centroid(cont_face, delt_nod)

switch cont_face.type
  case 'tria3' ,   center = [1/3,1/3,0]'; % 2D *******************************************
  case 'quad4' ,   center = [0,0,0]';
  otherwise   , error('Works only for linear 2D elements (tria3 or quad4).');
end
% Sizes
nele = cont_face.info.nele;
shapef = cont_face.fe.N(center);

delt_centroid = zeros(nele,3);
for i=1:nele
  delt_centroid(i,:) = shapef'*delt_nod(cont_face.nod(i,:),:);
end
