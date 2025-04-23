function test_linear_normals_centroids(cf_sl, cf_mast, delt_d_s, delt_d_m)
% tests linearizations of normals, tangents and centroids
% input argument is reference configuration

nN_s = cf_sl.info.nN; 
nele_s = cf_sl.info.nele;  
nN_ele_s = cf_sl.info.nN_ele;

nN_m = cf_mast.info.nN; 
nele_m = cf_mast.info.nele;  
nN_ele_m = cf_mast.info.nN_ele; 

xi_id_s = randi(nN_ele_s);
ele_id_s = randi(nele_s);
nod_id_s= randi(nN_s);

cf_sl_1 = cf_sl;  % cont_face.coo = X
cf_sl_2 = cf_sl;

cf_mast_1 = cf_mast;  % cont_face.coo = X
cf_mast_2 = cf_mast;

normals_storage_1 = get_normals_centroids(cf_sl_1); % current conf. normals centroids tangents

% derivatives from matrices -------
% ele normals
Nhatje = linear_get_nodal_element_normals(cf_sl_1, xi_id_s, ele_id_s); %Nhatje(x)
delt_nje_mat = Nhatje*reshape(delt_d_s',[3*nN_s, 1]);
ele_normal_1 = normals_storage_1.ele_normals{ele_id_s}(xi_id_s,:)';

% summed normals
Nhatj = linear_get_sum_normals(cf_sl_1, nod_id_s, normals_storage_1.ele_normals);
delt_nhatj_mat = Nhatj*reshape(delt_d_s',[3*nN_s, 1]);
sum_normal_1 = normals_storage_1.sum_normals(nod_id_s,:)';

% Averaged normals
Nj = linear_averaged_normals(cf_sl_1, normals_storage_1, nod_id_s);
delt_nj_mat = Nj*reshape(delt_d_s',[3*nN_s, 1]);
avg_normal_1 = normals_storage_1.averaged_normals(nod_id_s,:)';

% Tangents
[T1, T2] = linear_get_tangents(cf_sl_1, normals_storage_1, nod_id_s);
delt_t1_mat = T1*reshape(delt_d_s',[3*nN_s, 1]);
delt_t2_mat = T2*reshape(delt_d_s',[3*nN_s, 1]);
tangent1_1 = normals_storage_1.t1(nod_id_s,:)';
tangent2_1 = normals_storage_1.t2(nod_id_s,:)';

% n0
N0_s = linear_get_n0(cf_sl_1, normals_storage_1, ele_id_s);
delt_n0_mat = N0_s*reshape(delt_d_s',[3*nN_s, 1]);
normal0_1 = normals_storage_1.n0(ele_id_s,:)';

% x0
C0_s = linear_get_centroid(cf_sl_1, ele_id_s);
delt_x0_mat = C0_s*reshape(delt_d_s',[3*nN_s, 1]);
x0_1 = normals_storage_1.x0(ele_id_s,:)';

% add master
delt_d_c = [reshape(delt_d_m',[3*nN_m, 1]);reshape(delt_d_s',[3*nN_s, 1])];
N0 = [zeros(3,3*nN_m), N0_s]; C0 = [zeros(3,3*nN_m), C0_s];

% R
Ralg = linear_rotation_matrix(normal0_1, N0);
delt_R_mat = zeros(3);
for i=1:3
  for j=1:3
    delt_R_mat(i,j) = squeeze(Ralg(i,j,:))'*delt_d_c;
  end
end
sl.R_1 = rotation_matrix(normal0_1);


% init slave element
sl.nod_id = cf_sl.nod(ele_id_s,:);                      % row vec of glob nod indices in cf_sl
sl.coo_1 = cf_sl_1.coo(sl.nod_id, :)';                  % mat where each column are coo of nods in s
sl.x0_1 = x0_1;                                           % 3x1 coo of x0
sl.n0_1 = normal0_1;                                      % 3x1 coo of n0
sl.n = normals_storage_1.averaged_normals(sl.nod_id,:)';% mat where each column is avg normal
sl.fe = cf_sl.fe;                                       % sl finite element
% projections
sl.proj_1 = project_onto_plane(sl.coo_1, sl.x0_1, sl.n0_1);
% rotations
sl.rot_1 = (sl.R_1*(sl.proj_1 -sl.x0_1)) + sl.x0_1;

% dual shapef
[Ae_1, Me_1] = get_dual_shapef(sl.coo_1, sl.fe);
Aalg = linear_get_dual_shapef(sl.coo_1, sl.nod_id, sl.fe, Ae_1, Me_1, nN_s, nN_ele_s);
delt_A_mat = zeros(nN_ele_s);
for i=1:nN_ele_s
  for j=1:nN_ele_s
    delt_A_mat(i,j) = squeeze(Aalg(i,j,:))'*reshape(delt_d_s',[3*nN_s, 1]);
  end
end


for i=1:30
  ele_id_m = randi(nele_m);

  mast.nod_id = cf_mast.nod(ele_id_m,:);           % row vec of glob nod indices in cf_mast
  mast.coo_1 = cf_mast_1.coo(mast.nod_id, :)';  % mat where each column are coo of nods in m
  mast.fe = cf_mast.fe;                     % mast finite element
  mast.proj_1 = project_onto_plane(mast.coo_1, sl.x0_1, sl.n0_1);
  mast.rot_1 = (sl.R_1*(mast.proj_1 -sl.x0_1)) + sl.x0_1;
  
  % clipping in 2D with the vertex origins
  [rot_clip_1, clip_origin] = clipping_2D(sl.rot_1, mast.rot_1);
  if ~isempty(rot_clip_1)
    break;
  end
end

if isempty(rot_clip_1)
  fprint('No intersection found in 30 random attempts, try again');
  return;
end
rot_centr_1 = mean(rot_clip_1, 2);
% rotate polygon vertices back
clip_1 = (sl.R_1' * (rot_clip_1-sl.x0_1)) + sl.x0_1;
clip_centr_1 = (sl.R_1' * (rot_centr_1-sl.x0_1)) + sl.x0_1;

V = linear_rot_polyg_vert(sl.coo_1, mast.coo_1, sl.x0_1, sl.n0_1, N0, C0, ...
  clip_origin, sl.proj_1, mast.proj_1, sl.rot_1, mast.rot_1, nN_m, sl.R_1, Ralg, rot_clip_1);

ncells = size(V,1);
cell_id = randi(ncells);
delt_vert1_mat = V{cell_id,1}*delt_d_c;
delt_vert2_mat = V{cell_id,2}*delt_d_c;
delt_vert3_mat = V{cell_id,3}*delt_d_c;

vert1_1 = clip_centr_1;
vert2_1 = clip_1(:,cell_id);
vert3_1 = clip_1(:,mod(cell_id, ncells)+1);





% for greater error threshold the limits at epsilon>10^8 sometimes fail (prop floating point error)
tol = 10^-6;
for i=7:8
  epsilon = 10^-i;
  
  cf_sl_2.coo = cf_sl.coo + epsilon*delt_d_s;
  cf_mast_2.coo = cf_mast.coo + epsilon*delt_d_m;
  normals_storage_2 = get_normals_centroids(cf_sl_2);
  sl.coo_2 = cf_sl_2.coo(sl.nod_id, :)';
  mast.coo_2 = cf_mast_2.coo(mast.nod_id, :)'; 
  
% ---------- test ele normals ------------
  
  ele_normal_2 = normals_storage_2.ele_normals{ele_id_s}(xi_id_s,:)';
  delt_nje = (ele_normal_2-ele_normal_1)/epsilon;
  if (norm(delt_nje-delt_nje_mat) < tol)
     fprintf('For epsilon = %g: Ele normal at ele %d local node %d - Success\n', epsilon, ele_id_s, xi_id_s);
  else
     fprintf('For epsilon = %g: Ele normal at ele %d local node %d - Fail\n', epsilon, ele_id_s, xi_id_s);
  end
  
  
% ---------- test summed normals ----------
  
  sum_normal_2 = normals_storage_2.sum_normals(nod_id_s,:)';
  delt_nhatj = (sum_normal_2-sum_normal_1)/epsilon;
  if (norm(delt_nhatj-delt_nhatj_mat) < tol)
     fprintf('For epsilon = %g: Sum normal at node %d - Success\n', epsilon, nod_id_s);
  else
     fprintf('For epsilon = %g: Sum normal at node %d - Fail\n', epsilon, nod_id_s);
  end
  

% --------- test averaged normals ------
  
  avg_normal_2 = normals_storage_2.averaged_normals(nod_id_s,:)';
  delt_nj = (avg_normal_2-avg_normal_1)/epsilon;
  if (norm(delt_nj-delt_nj_mat) < tol)
     fprintf('For epsilon = %g: Avg normal at node %d - Success\n', epsilon, nod_id_s);
  else
     fprintf('For epsilon = %g: Avg normal at node %d - Fail\n', epsilon, nod_id_s);
  end
  
% --------- test tangents ------
  
  tangent1_2 = normals_storage_2.t1(nod_id_s,:)';
  delt_t1 = (tangent1_2-tangent1_1)/epsilon;
  if (norm(delt_t1-delt_t1_mat) < tol)
     fprintf('For epsilon = %g: Tangent 1 at node %d - Success\n', epsilon, nod_id_s);
  else
     fprintf('For epsilon = %g:  Tangent 1 at node %d - Fail\n', epsilon, nod_id_s);
  end

  tangent2_2 = normals_storage_2.t2(nod_id_s,:)';
  delt_t2 = (tangent2_2-tangent2_1)/epsilon;
  if (norm(delt_t2-delt_t2_mat) < tol)
     fprintf('For epsilon = %g: Tangent 2 at node %d - Success\n', epsilon, nod_id_s);
  else
     fprintf('For epsilon = %g:  Tangent 2 at node %d - Fail\n', epsilon, nod_id_s);
  end
  
% --------------- test n0 -----------------
  
  normal0_2 = normals_storage_2.n0(ele_id_s,:)';
  delt_n0 = (normal0_2-normal0_1)/epsilon;
  if (norm(delt_n0-delt_n0_mat) < tol)
     fprintf('For epsilon = %g: n0 normal at ele %d - Success\n', epsilon, ele_id_s);
  else
     fprintf('For epsilon = %g: n0 normal at ele %d - Fail\n', epsilon, ele_id_s);
  end
  
  % -------------- test x0 -----------------
  x0_2 = normals_storage_2.x0(ele_id_s,:)';
  delt_x0 = (x0_2-x0_1)/epsilon;
  if (norm(delt_x0-delt_x0_mat) < tol)
     fprintf('For epsilon = %g: x0 centroid at ele %d - Success\n', epsilon, ele_id_s);
  else
     fprintf('For epsilon = %g: x0 centroid at ele %d - Fail\n', epsilon, ele_id_s);
  end
  

%------------- test R -------------------
  R_2 = rotation_matrix(normal0_2);
  delt_R = (R_2-sl.R_1)/epsilon;
  if all(abs(delt_R(:) - delt_R_mat(:)) < tol)
    fprintf('For epsilon = %g: R at ele %d - Success\n', epsilon, ele_id_s);
  else
    fprintf('For epsilon = %g: R at ele %d - Fail\n', epsilon, ele_id_s);
  end
  
%------------- test dual shapef -------------------
  
  [Ae_2, Me_2] = get_dual_shapef(sl.coo_2, sl.fe);
  delt_A = (Ae_2-Ae_1)/epsilon;
  if all(abs(delt_A(:) - delt_A_mat(:)) < tol)
    fprintf('For epsilon = %g: A at ele %d - Success\n', epsilon, ele_id_s);
  else
    fprintf('For epsilon = %g: A at ele %d - Fail\n', epsilon, ele_id_s);
  end

% ---------- test clipping ----------------------
  
  % projections
  sl.proj_2 = project_onto_plane(sl.coo_2, x0_2, normal0_2);
  mast.proj_2 = project_onto_plane(mast.coo_2, x0_2, normal0_2);
  
  % rotations
  sl.rot_2 = (R_2*(sl.proj_2 -x0_2)) + x0_2;
  mast.rot_2 = (R_2*(mast.proj_2 -x0_2)) + x0_2;
  
  % clipping in 2D with the vertex origins
  [rot_clip_2, clip_origin] = clipping_2D(sl.rot_2, mast.rot_2);
  if isempty(rot_clip_2)
    fprintf('The clip was empty - no intersection of elements %d %d \n', ele_id_s, ele_id_m);
    return;
  end
  rot_centr_2 = mean(rot_clip_2, 2);
  % rotate polygon vertices back
  clip_2 = (R_2' * (rot_clip_2-x0_2)) + x0_2;
  clip_centr_2 = (R_2' * (rot_centr_2-x0_2)) + x0_2;

  vert1_2 = clip_centr_2;
  vert2_2 = clip_2(:,cell_id);
  vert3_2 = clip_2(:,mod(cell_id, ncells)+1);
  delt_vert1 = (vert1_2-vert1_1)/epsilon;
  delt_vert2 = (vert2_2-vert2_1)/epsilon;
  delt_vert3 = (vert3_2-vert3_1)/epsilon;
  if (norm(delt_vert1-delt_vert1_mat) < tol)
     fprintf('For epsilon = %g: v1 - centroid vertex at s ele %d and m ele %d and cell %d - Success\n', epsilon, ele_id_s, ele_id_m, cell_id);
  else
     fprintf('For epsilon = %g: v1 - centroid vertex at s ele %d and m ele %d and cell %d - Fail\n', epsilon, ele_id_s, ele_id_m, cell_id);
  end
  if (norm(delt_vert2-delt_vert2_mat) < tol)
     fprintf('For epsilon = %g: v2 at s ele %d and m ele %d and cell %d - Success\n', epsilon, ele_id_s, ele_id_m, cell_id);
  else
     fprintf('For epsilon = %g: v2 at s ele %d and m ele %d and cell %d - Fail\n', epsilon, ele_id_s, ele_id_m, cell_id);
  end
  if (norm(delt_vert3-delt_vert3_mat) < tol)
     fprintf('For epsilon = %g: v3 at s ele %d and m ele %d and cell %d - Success\n', epsilon, ele_id_s, ele_id_m, cell_id);
  else
     fprintf('For epsilon = %g: v3 at s ele %d and m ele %d and cell %d - Fail\n', epsilon, ele_id_s, ele_id_m, cell_id);
  end


end















