function [centroids, normals] = get_centroids_normals(mesh, face)
% given the mesh and face (as a number) computes normals of its
% elements 2D elements

  switch mesh.bou{face}.type
    case 'tria3' ,   center = [1/3,1/3,0]'; % 2D *******************************************
    case 'quad4' ,   center = [0,0,0]';
    otherwise   , error('Works only for linear 2D elements.');
  end

    nele = size(mesh.bou{face}.ele, 1); % Number of Elements on the face
    nN = size(mesh.bou{face}.nod, 2); % Number of Nodes on a single element

    fe = fe_init(mesh.bou{face}.type, 2); % Init finite element
    der_shapef = fe.dNdxi(center); % Deriv shape funcs of isoparametric (reference) element
          % we use center, but the derivatives are constant (if shape funcs are linear)
    shapef =  fe.N(center); % Shape funcs at centroid of isoparametric (reference) element

    t1 = zeros(3, nele); 
    t2 = zeros(3, nele);
    centroids = zeros(3, nele);
    
    for i=1:nele
      ele_nod_coo = mesh.nod.coo(mesh.bou{face}.nod(i,:), :); % coordinates of nodes of elements
        for j=1:nN
          t1(:,i) = t1(:,i) + der_shapef(j,1)*ele_nod_coo(j,:)'; % tangents in xi direction
          t2(:,i) = t2(:,i) + der_shapef(j,2)*ele_nod_coo(j,:)'; % tangents in eta direction
          centroids(:,i) = centroids(:,i) + shapef(j)*ele_nod_coo(j,:)';
        end
    end

    normt1 = sqrt(sum(t1.^2,1)); t1 = t1./normt1; % normed for each element
    normt2 = sqrt(sum(t2.^2,1)); t2 = t2./normt2;
    normals = cross(t2,t1,1); % needed to swap for right direction
end
