function Ctilde = linear_dual_mortar_fem(cont_face_s, cont_face_m, clips_storage, normals_storage)

% Sizes
nele_s = cont_face_s.info.nele;
nele_m = cont_face_m.info.nele;
nN_ele_s = cont_face_s.info.nN_ele;
nN_ele_m = cont_face_m.info.nN_ele;
nN_s = cont_face_s.info.nN;
nN_m = cont_face_m.info.nN; 

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
for i=1:nN_s
  averaged_normals(i,:) = sum_normals(i,:)/norm(sum_normals(i,:));
end
n0 = normals_storage.n0;
centroids = normals_storage.x0;

for k=1:nele_s
  s.nod = cont_face_s.nod(k,:);
  s.coo = cont_face_s.coo(s.nod, :);
  s.type = cont_face_s.type;
  s.x0 = centroids(k,:);
  s.n0 = n0(k,:);
  s.normals = averaged_normals(s.nod,:);
  s.fe = cont_face_s.fe;
  N0 = linear_get_n0(cont_face_s, normals_storage, k);    % get the 3xnN_s matrix for centroid normal
  C0 = linear_get_centroid(cont_face_s, k);               % get the 3xnN_s matrix for centroid
  N0 = [zeros(3,3*nN_m), N0]; C0 = [zeros(3,3*nN_m), C0]; % expand the matrices by a master part
  
  for i=1:nele_m
    if isempty(clips_storage{k,i})
      continue;
    end

    m.nod = cont_face_m.nod(i,:);
    m.coo = cont_face_m.coo(m.nod, :);
    m.fe = cont_face_m.fe;
    m.type = cont_face_m.type;

    clip = clips_storage{k,i}.vert; clip_origin = clips_storage{k,i}.orig;
    proj_s = clips_storage{k,i}.proj_s; proj_m = clips_storage{k,i}.proj_m;

    T_cell = linear_integr_cell(s.coo, m.coo, s.x0, s.n0, N0, C0, clip_origin, proj_s, proj_m, nN_m);

    Ae = clips_storage{k,i}.dshpf.Ae; Me = clips_storage{k,i}.dshpf.Me;
    delt_A = linear_get_dual_shapef(s.coo, s.delt_x, s.fe, Ae, Me);

    ncells = size(clip,1);
    clip_centr = mean(clip, 1);
    for j=1:ncells
      curr_cell = [clip_centr; clip(j,:); clip(mod(j, ncells)+1,:)];  % take a triangle segment
      delt_curr_cell = delt_integr_cells{j};
      delt_cell_gp = N_in_gp'*delt_curr_cell;    % deriv gp global coordinates on the segment (Popp 3D A16)
      s_proj_gp = clips_storage{k, i}.gp.proj_s{j}.coo;
      m_proj_gp = clips_storage{k, i}.gp.proj_m{j}.coo;
      s_alpha = clips_storage{k, i}.gp.proj_s{j}.alpha;
      m_alpha = clips_storage{k, i}.gp.proj_m{j}.alpha;
      s_delt_proj_gp = linear_proj_gp(s.coo, s.delt_x, s.fe, s_proj_gp, s_alpha, s.n0, s.delt_n0, delt_cell_gp); % [delt_xi, delt_eta, delt_alpha]
      m_delt_proj_gp = linear_proj_gp(m.coo, m.delt_x, m.fe, m_proj_gp, m_alpha, s.n0, s.delt_n0, delt_cell_gp);

      Ns_in_sgp = s.fe.N(s_proj_gp');
      Nm_in_mgp = m.fe.N(m_proj_gp');
      disc_gapf = dot(Ns_in_sgp'*s.normals, Nm_in_mgp'*m.coo-Ns_in_sgp'*s.coo, 2); % (A9) Popp 3D 
      delt_disc_gapf = linear_disc_gapf(s.coo, m.coo, s.delt_x, m.delt_x, s.fe, m.fe, s.delt_averaged_normals, ...
        s.normals, s_delt_proj_gp(:,1:2), m_delt_proj_gp(:,1:2), s_proj_gp, m_proj_gp);
      Jcell = clips_storage{k, i}.Jcell{j};
      delt_Jcell = linear_Jcell(curr_cell, delt_curr_cell);
      dPhidd = delt_A*Ns_in_sgp;   % dual shapef linear by displacement Popp 2D (32,33) Popp 3D (A22) 
      dshpf_in_sgp = Ae*Ns_in_sgp;
      
      
      for g=1:n_gp
        dNsdxi_in_sgp = s.fe.dNdxi(s_proj_gp(g,:)');
        dNmdxi_in_mgp = m.fe.dNdxi(m_proj_gp(g,:)');
        dPhidgp = Ae*dNsdxi_in_sgp; % dual shapef deriv by gauss point (nodal coo)

        delt_Phi = dPhidd(:,g) + dPhidgp*s_delt_proj_gp(g,1:2)';
        delt_weight_gap_tmp = delt_weight_gap_tmp + delt_Phi*disc_gapf(g)*w_gp(g)*Jcell + ...
          dshpf_in_sgp(:,g)*delt_disc_gapf(g)*w_gp(g)*Jcell + dshpf_in_sgp(:,g)*disc_gapf(g)*w_gp(g)*delt_Jcell;
        delt_D_tmp = delt_D_tmp + w_gp(g)*Jcell*dNsdxi_in_sgp*s_delt_proj_gp(g,1:2)'+Ns_in_sgp(:,g)*w_gp(g)*delt_Jcell;
        delt_M_tmp = delt_M_tmp + delt_Phi*Nm_in_mgp(:,g)'*w_gp(g)*Jcell + ...
          dshpf_in_sgp(:,g)*m_delt_proj_gp(g,1:2)*dNmdxi_in_mgp'*w_gp(g)*Jcell + ...
          dshpf_in_sgp(:,g)*Nm_in_mgp(:,g)'*w_gp(g)*delt_Jcell;
      end  
    end

    % map the local matrices into the "global" ones( "global" is still based on local nodal indexes on the face)
    delt_D(s.nod, s.nod) = delt_D(s.nod, s.nod) + diag(delt_D_tmp);
    delt_M(s.nod, m.nod) = delt_M(s.nod, m.nod) + delt_M_tmp;
    delt_weight_gap(s.nod) = delt_weight_gap(s.nod) + delt_weight_gap_tmp;
  end
end
mort_delt_D = nodal_blocks_to_mort_mat(delt_D);
mort_delt_M = nodal_blocks_to_mort_mat(delt_M);