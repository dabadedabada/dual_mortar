function test_linear_normals_centroids(cont_face)
% tests linearizations of normals, tangents and centroids
% input argument is reference configuration

nN = cont_face.info.nN; 
nele = cont_face.info.nele;  
nN_ele = cont_face.info.nN_ele; 

xi_id = randi(nN_ele);
ele_id = randi(nele);
nod_id = randi(nN);

err = 10^-6;
d = rand(nN,3);
delt_d = rand(nN,3);

cont_face_1 = cont_face;  % cont_face.coo = X
cont_face_2 = cont_face;

cont_face_1.coo = cont_face.coo + d; % cont_face_1.coo = x
normals_storage_1 = get_normals_centroids(cont_face_1); % current conf. normals centroids tangents

% derivatives from matrices -------
% ele normals
Nhatje = linear_get_nodal_element_normals(cont_face_1, xi_id, ele_id); %Nhatje(x)
delt_nje_mat = Nhatje*reshape(delt_d',[3*nN, 1]);

% summed normals
Nhatj = linear_get_sum_normals(cont_face_1, nod_id, normals_storage_1.ele_normals);
delt_nhatj_mat = Nhatj*reshape(delt_d',[3*nN, 1]);

% Averaged normals
Nj = linear_averaged_normals(cont_face_1, normals_storage_1, nod_id);
delt_nj_mat = Nj*reshape(delt_d',[3*nN, 1]);

% Tangents
[T1, T2] = linear_get_tangents(cont_face_1, normals_storage_1, nod_id);
delt_t1_mat = T1*reshape(delt_d',[3*nN, 1]);
delt_t2_mat = T2*reshape(delt_d',[3*nN, 1]);

% n0
N0 = linear_get_n0(cont_face_1, normals_storage_1, ele_id);
delt_n0_mat = N0*reshape(delt_d',[3*nN, 1]);

% x0
C0 = linear_get_centroid(cont_face_1, ele_id);
delt_x0_mat = C0*reshape(delt_d',[3*nN, 1]);

for i=4:8
  epsilon = 10^-i;
  
  cont_face_2.coo = cont_face.coo + d + epsilon*delt_d;
  normals_storage_2 = get_normals_centroids(cont_face_2);
  
  % ---------- test ele normals ------------
  %{
  n2 = normals_storage_2.ele_normals{ele_id}(xi_id,:)';
  n1 = normals_storage_1.ele_normals{ele_id}(xi_id,:)';
  delt_nje = (n2-n1)/epsilon;
  if (norm(delt_nje-delt_nje_mat) < err)
     fprintf('For epsilon = %g: Ele normal - Success\n', epsilon);
  else
     fprintf('For epsilon = %g: Ele normal - Fail\n', epsilon);
  end
  %}
  
  % ---------- test summed normals ----------
  %{
  n2 = normals_storage_2.sum_normals(nod_id,:)';
  n1 = normals_storage_1.sum_normals(nod_id,:)';
  delt_nhatj = (n2-n1)/epsilon;
  if (norm(delt_nhatj-delt_nhatj_mat) < err)
     fprintf('For epsilon = %g: Sum normal - Success\n', epsilon);
  else
     fprintf('For epsilon = %g: Sum normal - Fail\n', epsilon);
  end
  %}

  % --------- test averaged normals ------
  %{
  n2 = normals_storage_2.averaged_normals(nod_id,:)';
  n1 = normals_storage_1.averaged_normals(nod_id,:)';
  delt_nj = (n2-n1)/epsilon;
  if (norm(delt_nj-delt_nj_mat) < err)
     fprintf('For epsilon = %g: Avg normal - Success\n', epsilon);
  else
     fprintf('For epsilon = %g: Avg normal - Fail\n', epsilon);
  end
  %}
  % --------- test tangents ------
  %{
  t12 = normals_storage_2.t1(nod_id,:)';
  t11 = normals_storage_1.t1(nod_id,:)';
  delt_t1 = (t12-t11)/epsilon;
  if (norm(delt_t1-delt_t1_mat) < err)
     fprintf('For epsilon = %g: Tangent 1 - Success\n', epsilon);
  else
     fprintf('For epsilon = %g:  Tangent 1 - Fail\n', epsilon);
  end

  t22 = normals_storage_2.t2(nod_id,:)';
  t21 = normals_storage_1.t2(nod_id,:)';
  delt_t2 = (t22-t21)/epsilon;
  if (norm(delt_t2-delt_t2_mat) < err)
     fprintf('For epsilon = %g: Tangent 2 - Success\n', epsilon);
  else
     fprintf('For epsilon = %g:  Tangent 2 - Fail\n', epsilon);
  end
  %}
  % --------------- test n0 -----------------
  
  n2 = normals_storage_2.n0(ele_id,:)';
  n1 = normals_storage_1.n0(ele_id,:)';
  delt_n0 = (n2-n1)/epsilon;
  if (norm(delt_n0-delt_n0_mat) < err)
     fprintf('For epsilon = %g: n0 normal - Success\n', epsilon);
  else
     fprintf('For epsilon = %g: n0 normal - Fail\n', epsilon);
  end
  
  % -------------- test x0 -----------------
  x02 = normals_storage_2.x0(ele_id,:)';
  x01 = normals_storage_1.x0(ele_id,:)';
  delt_x0 = (x02-x01)/epsilon;
  if (norm(delt_x0-delt_x0_mat) < err)
     fprintf('For epsilon = %g: x0 normal - Success\n', epsilon);
  else
     fprintf('For epsilon = %g: x0 normal - Fail\n', epsilon);
  end
end















