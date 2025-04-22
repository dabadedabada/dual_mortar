function Ctilde = linear_dual_mortar_fem(cf_sl, cf_mast, clips_storage, normals_storage, slave_storage, master_storage)

% Sizes
nele_s = cf_sl.info.nele;
nele_m = cf_mast.info.nele;
nN_ele_s = cf_sl.info.nN_ele;
nN_ele_m = cf_mast.info.nN_ele;
nN_s = cf_sl.info.nN;
nN_m = cf_mast.info.nN; 

% For quad element or distorted geometry choose high order,
% here o=7 is pretty accurate (13 points)
fe_cell = fe_init('tria3', 7);   % Init finite element for cells
gauss_points = fe_cell.QP;       % gauss points
w_gp = fe_cell.QW;               % weights
n_gp = size(gauss_points, 2);     % number of gp
N_in_gp = fe_cell.N(gauss_points);  % shape funcs in gauss points

% init linearizations for D,M,weight_gap
delt_D = zeros(nN_s);
delt_M = zeros(nN_s, nN_m);
delt_weight_gap = zeros(nN_s,1);

sum_normals = normals_storage.sum_normals;
averaged_normals = zeros(size(sum_normals));
for m=1:nN_s
  averaged_normals(m,:) = sum_normals(m,:)/norm(sum_normals(m,:));
end

for s=1:nele_s
  sl.nod = cf_sl.nod(s,:);
  sl.coo = cf_sl.coo(sl.nod, :)';
  sl.type = cf_sl.type;
  sl.x0 = normals_storage.x0(s,:)';
  sl.n0 = normals_storage.n0(s,:)';
  sl.normals = averaged_normals(sl.nod,:)';
  sl.fe = cf_sl.fe;
  N0 = linear_get_n0(cf_sl, normals_storage, s);    % get the 3xnN_s matrix for centroid normal
  C0 = linear_get_centroid(cf_sl, s);               % get the 3xnN_s matrix for centroid
  N0 = [zeros(3,3*nN_m), N0]; C0 = [zeros(3,3*nN_m), C0]; % expand the matrices by a master part

  sl.R = slave_storage{s}.R;
  Ralg = linear_rotation_matrix(sl.n0, N0);
  sl.proj = slave_storage{s}.proj;
  sl.rot = slave_storage{s}.rot;

  Ae = slave_storage{s}.Ae; 
  Me = slave_storage{s}.Me;

  Aalg = linear_get_dual_shapef(sl.coo, sl.nod, sl.fe, Ae, Me, nN_s, nN_ele_s);
  
  for m=1:nele_m
    if isempty(clips_storage{s,m})
      continue;
    end

    mast.nod = cf_mast.nod(m,:);
    mast.coo = cf_mast.coo(mast.nod, :)';
    mast.fe = cf_mast.fe;
    mast.type = cf_mast.type;

    mast.proj = master_storage{s}.proj;
    mast.rot = master_storage{s}.rot;

    rot_clip = clips_storage{s,m}.rot_vert;
    clip = clips_storage{s,m}.vert; 
    clip_origin = clips_storage{s,m}.orig;
    rot_centr = mean(rot_clip, 2);
    clip_centr = (sl.R' * (rot_centr-sl.x0)) + sl.x0;
    
    % linearization of projection to aux plane (P), rotation (Ptilde), 2D clipping (Vtilde) and
    % roation back (V) - ncells x 3 cell of matrices V of cell vertices
    V = linear_rot_polyg_vert(sl.coo, mast.coo, sl.x0, sl.n0, N0, C0, ...
      clip_origin, sl.proj, mast.proj, sl.rot, mast.rot, nN_m, sl.R, Ralg, rot_clip);

    Dalg_loc = zeros(nN_ele_s,nN_ele_s,3*(nN_m + nN_s));
    Malg_loc = zeros(nN_ele_s,nN_ele_m,3*(nN_m + nN_s));


    ncells = size(clip,1);
    for cell_id=1:ncells
      % take a triangle segment
      cell_vert_coo = [clip_centr, clip(:,cell_id), clip(:,mod(cell_id, ncells)+1)]; % 3x3 [v1,v2,v3]
      % gp global coordinates on the segment (Popp diss A.30)
      curr_cell_gp = cell_vert_coo*N_in_gp;  % 3xn_gp
      J_cell = clips_storage{s, m}.Jcell{cell_id};
      jcell = linear_cell_Jacobian(cell_vert_coo, V{cell_id,1}, V{cell_id,2}, V{cell_id,3});
      for gp=1:n_gp
        Ghat = linear_cell_gp(V{cell_id,1}, V{cell_id,2}, V{cell_id,3}, fe_cell, gp);
        s_proj_gp = clips_storage{s, m}.gp.sl.proj{cell_id}.coo(:,gp);   % sl   [xi; eta]
        alpha = clips_storage{s, m}.gp.sl.proj{cell_id}.alpha(gp);     % alpha
        m_proj_gp = clips_storage{s, m}.gp.mast.proj{cell_id}.coo(:,gp); % mast [xi; eta]
        beta = clips_storage{s, m}.gp.mast.proj{cell_id}.alpha(gp);    % beta (alpha for master)
        G_s = linear_proj_gp(sl.coo, sl.nod, sl.fe, s_proj_gp, alpha, sl.n0, N0, Ghat, 1, nN_m);
        G_m = linear_proj_gp(mast.coo, mast.nod, mast.fe, m_proj_gp, beta, sl.n0,N0, Ghat, 2, nN_m);
        F = linear_dshpf_in_sl_proj_gp(Ae, Aalg, sl.fe, s_proj_gp, G_s, nN_m);

        

      end  
    end
  end
end
