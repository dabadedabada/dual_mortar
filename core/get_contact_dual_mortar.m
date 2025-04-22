function [mort_D, mort_M, weight_gap, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cf_sl, cf_mast, normals_storage)

% Sizes
nele_s = cf_sl.info.nele;
nele_m = cf_mast.info.nele;
nN_ele_s = cf_sl.info.nN_ele;
nN_ele_m = cf_mast.info.nN_ele;
nN_s = cf_sl.info.nN;
nN_m = cf_mast.info.nN; 

% ----- choose number of integration points on integration cell -------
% For quad element or distorted geometry choose high order,
% here o=7 is pretty accurate (13 points)
fe_cell = fe_init('tria3', 7);   % Init finite element for cells
gauss_points = fe_cell.QP;       % gauss points
w_gp = fe_cell.QW;               % weights
n_gp = size(gauss_points, 2);     % number of gp
N_in_gp = fe_cell.N(gauss_points);  % shape funcs in gauss points

% init D,M before expanding and weight_gap in slave cont nodes
D = zeros(nN_s);
M = zeros(nN_s, nN_m);
weight_gap = zeros(nN_s,1);

% Averaged normals are normed sums in nodes, so we only hold the sums
sum_normals = normals_storage.sum_normals;
averaged_normals = zeros(size(sum_normals));
for m=1:nN_s
  averaged_normals(m,:) = sum_normals(m,:)/norm(sum_normals(m,:));
end

% init storage of clip polygons
clips_storage = cell(nele_s, nele_m);
slave_storage = cell(nele_s);
master_storage = cell(nele_m);

% Mortar coupling algorithm 
for s=1:nele_s
  % Init slave element
  sl.nod_id = cf_sl.nod(s,:);           % row vec of glob nod indices in cf_sl
  sl.coo = cf_sl.coo(sl.nod_id, :)';    % mat where each column are coo of nods in s
  sl.x0 = normals_storage.x0(s,:)';     % 3x1 coo of x0
  sl.n0 = normals_storage.n0(s,:)';     % 3x1 coo of n0
  sl.n = averaged_normals(sl.nod_id,:)';% mat where each column is avg normal
  sl.fe = cf_sl.fe;                     % sl finite element

  % project slave nodes onto aux plan - sl.proj
  sl.proj = project_onto_plane(sl.coo, sl.x0, sl.n0);

  % Rodriguez's rotation matrix for rot to horizontal plane  
  sl.R = rotation_matrix(sl.n0);
  
  % rotate slave points
  sl.rot = (sl.R*(sl.proj -sl.x0)) + sl.x0;
  
  % Dual shape functions matrices (Me for linearization later)
  [Ae, Me] = get_dual_shapef(sl.coo, sl.fe);

  slave_storage{s}.R = sl.R;
  slave_storage{s}.proj =sl.proj;
  slave_storage{s}.rot = sl.rot;
  slave_storage{s}.Ae = Ae;
  slave_storage{s}.Me = Me;

  for m=1:nele_m
    % Init master element
    mast.nod_id = cf_mast.nod(m,:);           % row vec of glob nod indices in cf_mast
    mast.coo = cf_mast.coo(mast.nod_id, :)';  % mat where each column are coo of nods in m
    mast.fe = cf_mast.fe;                     % mast finite element

    D_tmp = zeros(nN_ele_s);
    M_tmp = zeros(nN_ele_s,nN_ele_m);
    weight_gap_tmp = zeros(nN_ele_s,1);
    
    % project mast nodes onto aux plane and rotate to horizontal plane
    mast.proj = project_onto_plane(mast.coo, sl.x0, sl.n0);
    mast.rot = (sl.R*(mast.proj -sl.x0)) + sl.x0;

    master_storage{m}.proj = mast.proj;
    master_storage{m}.rot = mast.rot;

    %plot_3d_polygon(sl.rot', mast.rot');

    % clipping in 2D with the vertex origins
    [rot_clip, clip_origin] = clipping_2D(sl.proj, mast.proj);
    if isempty(rot_clip)
      continue; % if theres no intersection of elements, skip iteration
    end
    rot_centr = mean(rot_clip, 2);
    % rotate polygon vertices back
    clip = (sl.R' * (rot_clip-sl.x0)) + sl.x0;
    clip_centr = (sl.R' * (rot_centr-sl.x0)) + sl.x0;

    % Assign to storage for current slave and master element
    clips_storage{s, m}.rot_vert = rot_clip;
    clips_storage{s, m}.vert = clip;
    clips_storage{s, m}.orig = clip_origin;
    clips_storage{s, m}.sl.proj = sl.proj;
    clips_storage{s, m}.mast.proj = mast.proj;

    % average for clip centroid
    
    % ncells is equal to number of clip vertices
    ncells = size(clip, 1);
 
    for cell_id=1:ncells
      % take a triangle segment
      cell_vert_coo = [clip_centr, clip(:,cell_id), clip(:,mod(cell_id, ncells)+1)]; % 3x3 [v1,v2,v3]
      % gp global coordinates on the segment (Popp diss A.30)
      curr_cell_gp = cell_vert_coo*N_in_gp;  % 3xn_gp  
    
      %plot_points_3d(cell_gp');
      % project cell gp along n0 to sl and mast surface and find
      % isoparametric ele coordinates
      s_proj_gp = newton_it_for_projgp(sl.coo, sl.fe, sl.n0, curr_cell_gp);    % sl   3xn_gp [xi; eta; alpha]
      m_proj_gp = newton_it_for_projgp(mast.coo, mast.fe, sl.n0, curr_cell_gp);% mast 3xn_gp [xi; eta; beta]

      clips_storage{s, m}.gp.sl.proj{cell_id}.coo = s_proj_gp(1:2,:);   % sl   [xi; eta]
      clips_storage{s, m}.gp.sl.proj{cell_id}.alpha = s_proj_gp(3,:);   % alpha
      clips_storage{s, m}.gp.mast.proj{cell_id}.coo = m_proj_gp(1:2,:); % mast [xi; eta]
      clips_storage{s, m}.gp.mast.proj{cell_id}.alpha = m_proj_gp(3,:); % beta
    
      %plot_points_on_ref_element(s.type, s_proj_gp{i})
    
      % cell Jacobian (Popp diss A.24)
      J_cell = norm(cross(cell_vert_coo(:,2)-cell_vert_coo(:,1),cell_vert_coo(:,3)-cell_vert_coo(:,1)));
      clips_storage{s, m}.Jcell{cell_id} = J_cell;
    
      Ns_in_sgp = sl.fe.N(s_proj_gp(1:2,:));   % values of sl shpf at sl proj of gp: nN_ele_s x n_gp
      Nm_in_mgp = mast.fe.N(m_proj_gp(1:2,:)); % values of mast shpf at mast proj of gp: nN_ele_m x n_gp
    
      dshpf_in_sgp = Ae*Ns_in_sgp;             % values of dshpf at sl proj of gp: nN_ele_s x n_gp
       
      % outward normals at gp on curr conf sl surf
      normals_gp = normalize(sl.n*Ns_in_sgp,1);
      
      % discrete gaps between curr coo of sl and mast gp proj (A9) Popp 3D
      disc_gapf = dot(normals_gp, mast.coo*Nm_in_mgp-sl.coo*Ns_in_sgp, 1); % 1 x n_gp 
      
      % Mortar integrals (popp 4.35)
      % local D and M for one element pair (s,m)
      for gp=1:n_gp
        D_tmp = D_tmp + dshpf_in_sgp(:,gp)*Ns_in_sgp(:,gp)'*w_gp(gp)*J_cell; 
        M_tmp = M_tmp + dshpf_in_sgp(:,gp)*Nm_in_mgp(:,gp)'*w_gp(gp)*J_cell;
        weight_gap_tmp = weight_gap_tmp + dshpf_in_sgp(:,gp)*disc_gapf(gp)*w_gp(gp)*J_cell;
      end
    end
    
    % map the local matrices into the "global" ones( "global" is still based on local nodal indexes on the face)
    D(sl.nod_id, sl.nod_id) = D(sl.nod_id, sl.nod_id) + D_tmp;
    M(sl.nod_id, mast.nod_id) = M(sl.nod_id, mast.nod_id) + M_tmp;
    weight_gap(sl.nod_id) = weight_gap(sl.nod_id) + weight_gap_tmp;
  end
end

% here mort_D^(3n x 3n) and mort_M^(3n x 3m)
% Popp 3D (26) and (27)
mort_D = nodal_blocks_to_mort_mat(D);
mort_M = nodal_blocks_to_mort_mat(M);


