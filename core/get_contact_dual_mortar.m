function [mort_D, mort_M, weight_gap, clips_storage] = get_contact_dual_mortar(cont_face_s, cont_face_m, normals_storage)

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

% init D,M,weight_gap
D = zeros(nN_s);
M = zeros(nN_s, nN_m);
weight_gap = zeros(nN_s,1);

% Averaged normals are normed sums in nodes, so we only hold the sums
% Need both for linearizations
sum_normals = normals_storage.sum_normals;
averaged_normals = zeros(size(sum_normals));
for i=1:nN_s
  averaged_normals(i,:) = sum_normals(i,:)/norm(sum_normals(i,:));
end
n0 = normals_storage.n0;
centroids = normals_storage.x0;

% Algorithm 1 (Popp 3d) 
clips_storage = cell(nele_s, nele_m);

for k=1:nele_s
  % Init slave element
  s.nod = cont_face_s.nod(k,:);
  s.coo = cont_face_s.coo(s.nod, :);
  s.x0 = centroids(k,:);
  s.n0 = n0(k,:);
  s.normals = averaged_normals(s.nod,:);
  s.fe = cont_face_s.fe;

  proj_s = project_onto_plane(s.coo, s.x0, s.n0); % project onto aux plane

  for i=1:nele_m
    % Init master element
    m.nod = cont_face_m.nod(i,:);
    m.coo = cont_face_m.coo(m.nod, :);
    m.fe = cont_face_m.fe;

    D_tmp = zeros(nN_ele_s);
    M_tmp = zeros(nN_ele_s,nN_ele_m);
    weight_gap_tmp = zeros(nN_ele_s,1);
    
    proj_m = project_onto_plane(m.coo, s.x0, s.n0);

    % Rotate to 2D, clip and rotate back
    [clip, clip_origin] = aux_project_and_clip(s.n0, s.x0, proj_s, proj_m);
    if isempty(clip)
      continue; % if theres no intersection of elements, skip iteration
    end

    % Assign to storage for current slave and master element
    clips_storage{k, i}.vert = clip;
    clips_storage{k, i}.orig = clip_origin;
    clips_storage{k, i}.proj_s = proj_s;
    clips_storage{k, i}.proj_m = proj_m;

    % average for center
    clip_centr = mean(clip, 1);
    ncells = size(clip, 1);
 
    for j=1:ncells
      curr_cell = [clip_centr; clip(j,:); clip(mod(j, ncells)+1,:)];  % take a triangle segment
      cell_gp = N_in_gp'*curr_cell;    % gp global coordinates on the segment (Popp diss A.30)
    
      % project to slave and master isoparametric element along slave normal n0
      s_proj_gp = newton_it_for_projgp(s.coo, s.fe, s.n0, cell_gp);  
      m_proj_gp = newton_it_for_projgp(m.coo, m.fe, s.n0, cell_gp);

      clips_storage{k, i}.gp.global{j} = cell_gp;
      clips_storage{k, i}.gp.proj_s{j}.coo = s_proj_gp(:,1:2);
      clips_storage{k, i}.gp.proj_s{j}.alpha = s_proj_gp(:,3);
      clips_storage{k, i}.gp.proj_m{j}.coo = m_proj_gp(:,1:2);
      clips_storage{k, i}.gp.proj_m{j}.alpha = m_proj_gp(:,3);
    
      %plot_points_on_ref_element(s.type, s_proj_gp{i})
    
      % cell Jacobian (Popp diss A.24)
      J_cell = norm(cross(curr_cell(2,:)-curr_cell(1,:),curr_cell(3,:)-curr_cell(1,:), 2));
      clips_storage{k, i}.Jcell{j} = J_cell;

      % Dual shape functions matrices (Me for linearization later)
      [Ae, Me] = get_dual_shapef(s.coo,s.fe);
      clips_storage{k, i}.dshpf.Ae{j} = Ae;
      clips_storage{k, i}.dshpf.Me{j} = Me;
    
      Ns_in_sgp = s.fe.N(s_proj_gp(:,1:2)');
      Nm_in_mgp = m.fe.N(m_proj_gp(:,1:2)');
    
      dshpf_in_sgp = Ae*Ns_in_sgp;
      disc_gapf = dot(Ns_in_sgp'*s.normals, Nm_in_mgp'*m.coo-Ns_in_sgp'*s.coo, 2); % (A9) Popp 3D 
      % Ns_in_sgp'*s.normals should be normed, as its the n_gp normal ???
      
      % Mortar integrals (popp 4.35)
      % local D and M for one element pair (s,m)
      for g=1:n_gp
        D_tmp = D_tmp + dshpf_in_sgp(:,g)*Ns_in_sgp(:,g)'*w_gp(g)*J_cell; 
        M_tmp = M_tmp + dshpf_in_sgp(:,g)*Nm_in_mgp(:,g)'*w_gp(g)*J_cell;
        weight_gap_tmp = weight_gap_tmp + dshpf_in_sgp(:,g)*disc_gapf(g)*w_gp(g)*J_cell;
      end
    end
    
    % map the local matrices into the "global" ones( "global" is still based on local nodal indexes on the face)
    D(s.nod, s.nod) = D(s.nod, s.nod) + D_tmp;
    M(s.nod, m.nod) = M(s.nod, m.nod) + M_tmp;
    weight_gap(s.nod) = weight_gap(s.nod) + weight_gap_tmp;
  end
end

% here mort_D^(3n x 3n) and mort_M^(3n x 3m)
% Popp 3D (26) and (27)
mort_D = nodal_blocks_to_mort_mat(D);
mort_M = nodal_blocks_to_mort_mat(M);


