function mesh = generate_cylinder_mesh(r, h, ne, etype)


Ly = h;
Lx = r*sqrt(2);
Lz = Lx;
mesh = mesh_generator([Lx,Ly,Lz], ne, etype);

% number of nodes in any direction, also number of slices in that direction
nNx = ne(1) + 1;
nNy = ne(2) + 1;
nNz = ne(3) + 1;

% number of nodes on given slice
nN_slicex = nNy*nNz;
nN_slicey = nNx*nNz;
nN_slicez = nNx*nNy;

% make each slice in y-direction into a circle
%r = norm(mesh.nod.coo(1,:)-mesh.nod.coo(nNx+(nNz-1)*nN_slicey,:));
for j=1:nNy
  temp_n = nNz -2;
  temp_slicey_ids = zeros(1,2*temp_n);
  id1_slicey = 1 + (j-1)*nNx;
  for i=1:temp_n
    id = 1 + (i-1)*2;
    temp_slicey_ids(id) = id1_slicey + nN_slicez*i;
    temp_slicey_ids(id+1) = id1_slicey + nN_slicez*i + nNx-1;
  end
  idend_slicey = id1_slicey + nNx - 1 + (nNz-1)*nN_slicez;
  slicey_ids = [id1_slicey:id1_slicey+ nNx - 1, temp_slicey_ids, idend_slicey - (nNx - 1):idend_slicey];
  n = length(slicey_ids);
  for i=1:n/2
    A = mesh.nod.coo(slicey_ids(i),:);
    B = mesh.nod.coo(slicey_ids(n-(i-1)),:);
    [newpoint1, newpoint2] = extend_points_to_r(A,B,r);
    mesh.nod.coo(slicey_ids(i),:) = newpoint1;
    mesh.nod.coo(slicey_ids(n-(i-1)),:) = newpoint2;
  end
  %plot_points_3d(mesh.nod.coo(slicey_ids,:));
end

