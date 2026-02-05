clear
%addpath(genpath('..'));
addpath("core");
addpath("plot");
addpath("core_matfem");
addpath(genpath("linearizations"));
addpath("clipper2");

mesh  = cell(2,1);
etype = cell(2,1);
etype{1} = 'hexa8';
etype{2} = 'hexa8'; %'tetr4';

% Create reference configurations X^(1,2), t=0
mesh{1} = mesh_generator([3,3,2], [6,1,1],etype{1});
mesh{2} = mesh_generator([12,12,1], [2,2,1],etype{2});

%mesh{1}.nod.coo=mesh{1}.nod.coo + [0,0,1.1];
%mesh{1}.nod.coo=mesh{1}.nod.coo + [0.5,1.5,0];

mesh{1}.nod.coo = mesh{1}.nod.coo + [-1.5,-1.5,0 ];
mesh{2}.nod.coo = mesh{2}.nod.coo + [-6.0,-6.0,0 ];
if false
  mesh{1}.nod.coo = mesh{1}.nod.coo + [0,0,1.1];
else
  tmp = 4;
  d  = tmp - mesh{1}.nod.coo(:,3);
  al = atan2(mesh{1}.nod.coo(:,1),tmp-(0*mesh{1}.nod.coo(:,3)+min(mesh{1}.nod.coo(:,3))));
  mesh{1}.nod.coo = [d.*sin(al)-1 mesh{1}.nod.coo(:,2)  tmp-d.*cos(al)+1.01];
end

%plot_meshes(mesh);

% Init contact faces
cont_face = cell(2,1);
cont_face{1} = init_cont_face(mesh{1}, 5);
cont_face{2} = init_cont_face(mesh{2}, 6);

nN_s = cont_face{1}.info.nN;
nN_m = cont_face{2}.info.nN;

% Normals, centroids and tangents
normals_storage = get_normals_centroids(cont_face{1});

% --------------------------------------------------------------------------%
% ---------------- testing linearizations ----------------------------------%
epsilon = 10^-6;

% ----------- element normal directions --------------------------%
%{
for s=1:cont_face{1}.info.nele
  Dn_alg = zeros(3*cont_face{1}.info.nN_ele, 3*nN_s);
  Dn_lim = zeros(3*cont_face{1}.info.nN_ele, 3*nN_s);
  for xi_id=1:cont_face{1}.info.nN_ele
    Dn_alg(xi_id*3-2:xi_id*3,:) = linear_get_nodal_element_normals(cont_face{1}, xi_id, s);
  end
  for i=1:nN_s
    for d=1:3
      cont_face_2 = cont_face{1};
      cont_face_2.coo(i,d) = cont_face_2.coo(i,d) + epsilon;
      normals_storage_2 = get_normals_centroids(cont_face_2);
      Dn_lim(:,d + 3*(i-1)) = (reshape(normals_storage_2.ele_normals{s}', 3*cont_face{1}.info.nN_ele, 1)...
        - reshape(normals_storage.ele_normals{s}', 3*cont_face{1}.info.nN_ele, 1))/epsilon;
    end
  end
  fprintf('Element normal directions: %g\n',max(reshape(abs(Dn_lim-Dn_alg),[],1)));
end
%}

% --------------- tangents -------------------------------%
Dt1_alg = zeros(3*nN_s, 3*nN_s); Dt2_alg = zeros(3*nN_s, 3*nN_s);
Dt1_lim = zeros(size(Dt1_alg)); Dt2_lim = zeros(size(Dt1_alg));

for i=1:nN_s
  [T1, T2] = linear_get_tangents(cont_face{1}, normals_storage, i);
  Dt1_alg(3*i-2:3*i,:) = T1; Dt2_alg(3*i-2:3*i,:) = T2; 
  for d=1:3
    cont_face_2 = cont_face{1};
    cont_face_2.coo(i,d) = cont_face_2.coo(i,d) + epsilon;
    normals_storage_2 = get_normals_centroids(cont_face_2);
    Dt1_lim(:,d + 3*(i-1)) = (reshape(normals_storage_2.t1', 3*nN_s, 1)-reshape(normals_storage.t1', 3*nN_s, 1))/epsilon;
    Dt2_lim(:,d + 3*(i-1)) = (reshape(normals_storage_2.t2', 3*nN_s, 1)-reshape(normals_storage.t2', 3*nN_s, 1))/epsilon;
  end
end
fprintf('t1: %g\n',max(reshape(abs(Dt1_alg-Dt1_lim),[],1)));
fprintf('t2: %g\n',max(reshape(abs(Dt2_alg-Dt2_lim),[],1)));

%{
% -------------- element centroids, centroid normals ------------------------------%
Dn_alg = zeros(3*cont_face{1}.info.nele, 3*nN_s); Dn_lim = zeros(size(Dn_alg));
Dx0_alg = zeros(3*cont_face{1}.info.nele, 3*nN_s); Dx0_lim = zeros(size(Dx0_alg));

for s=1:cont_face{1}.info.nele
  Dn_alg(s*3-2:s*3,:) = linear_get_n0(cont_face{1}, normals_storage, s);
  Dx0_alg(s*3-2:s*3,:) = linear_get_centroid(cont_face{1}, s);
  for i=1:nN_s
    for d=1:3
      cont_face_2 = cont_face{1};
      cont_face_2.coo(i,d) = cont_face_2.coo(i,d) + epsilon;
      normals_storage_2 = get_normals_centroids(cont_face_2);
      Dn_lim(3*s-2:3*s,d + 3*(i-1)) = (normals_storage_2.n0(s,:)'-normals_storage.n0(s,:)')/epsilon;
      Dx0_lim(3*s-2:3*s,d + 3*(i-1)) = (normals_storage_2.x0(s,:)'-normals_storage.x0(s,:)')/epsilon;
    end
  end
end
fprintf('Centroids: %g\n',max(reshape(abs(Dx0_lim-Dx0_alg),[],1)));
fprintf('Centroid normals: %g\n',max(reshape(abs(Dn_lim-Dn_alg),[],1)));
%}
%{
% --------------- rotation matrix ------------------------------------%
DR_alg = zeros(cont_face{1}.info.nele, 3, 3, 3*nN_s);
DR_lim = zeros(size(DR_alg));
for s=1:cont_face{1}.info.nele
  DR_alg(s,:,:,:) = linear_rotation_matrix(normals_storage.n0(s,:)', linear_get_n0(cont_face{1}, normals_storage, s));
  for i=1:nN_s
    for d=1:3
      cont_face_2 = cont_face{1};
      cont_face_2.coo(i,d) = cont_face_2.coo(i,d) + epsilon;
      normals_storage_2 = get_normals_centroids(cont_face_2);
      DR_lim(s,:,:,d + 3*(i-1)) = (rotation_matrix(normals_storage_2.n0(s,:)')-rotation_matrix(normals_storage.n0(s,:)'))/epsilon;
    end
  end
end
fprintf('Rotation matrix: %g \n', max(reshape(abs(DR_lim-DR_alg),[],1)));
%}
%{
% ------------------ dual shpf ----------------------------------%
[mort_D, mort_M, weight_gap, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cont_face{1}, cont_face{2}, normals_storage);
DA_alg = zeros(cont_face{1}.info.nele, cont_face{1}.info.nN_ele, cont_face{1}.info.nN_ele, 3*nN_s);
DA_lim = zeros(size(DA_alg));
for s=1:cont_face{1}.info.nele
  DA_alg(s,:,:,:) = linear_get_dual_shapef(cont_face{1}.coo(cont_face{1}.nod(s,:),:)', cont_face{1}.nod(s,:), cont_face{1}.fe,...
    slave_storage{s}.Ae, slave_storage{s}.Me, nN_s, cont_face{1}.info.nN_ele);
  for i=1:nN_s
    for d=1:3
      cont_face_2 = cont_face{1};
      cont_face_2.coo(i,d) = cont_face_2.coo(i,d) + epsilon;
      normals_storage_2 = get_normals_centroids(cont_face_2);
      [mort_D, mort_M, weight_gap, clips_storage, slave_storage_2, master_storage] = get_contact_dual_mortar(cont_face_2, cont_face{2}, normals_storage_2);
      DA_lim(s,:,:,d + 3*(i-1)) = (slave_storage_2{s}.Ae-slave_storage{s}.Ae)/epsilon;
    end
  end
end
fprintf('Dshpf Ae matrix: %g \n', max(reshape(abs(DA_lim-DA_alg),[],1)));
%}
%{
% ------------------ projections and rotations to 2D -------------------------------%
[mort_D, mort_M, weight_gap, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cont_face{1}, cont_face{2}, normals_storage);
Dproj_alg = zeros(cont_face{1}.info.nele,3*cont_face{1}.info.nN_ele, 3*(nN_m+nN_s));
Dproj_lim = zeros(size(Dproj_alg));
Drot_alg = zeros(size(Dproj_alg)); Drot_lim = zeros(size(Dproj_lim));
Drotback_alg = zeros(size(Dproj_alg)); Drotback_lim = zeros(size(Dproj_lim));
Dclip_alg = zeros(size(Dproj_alg)); Dclip_lim = zeros(size(Dproj_lim));
for s=1:cont_face{1}.info.nele
  for xi_id=1:cont_face{1}.info.nN_ele
    Dproj_alg(s,xi_id*3-2:xi_id*3,:) = linear_project_onto_plane(cont_face{1}.coo(cont_face{1}.nod(s,xi_id),:)', cont_face{1}.nod(s,xi_id), normals_storage.x0(s,:)', normals_storage.n0(s,:)',...
      [zeros(3,3*nN_m),linear_get_n0(cont_face{1}, normals_storage, s)], [zeros(3,3*nN_m), linear_get_centroid(cont_face{1}, s)], 1, nN_m);
    Drot_alg(s,xi_id*3-2:xi_id*3,:) = linear_rotate_points(slave_storage{s}.R, linear_rotation_matrix(normals_storage.n0(s,:)', [zeros(3,3*nN_m),linear_get_n0(cont_face{1}, normals_storage, s)]),...
      slave_storage{s}.proj(:,xi_id), squeeze(Dproj_alg(s,xi_id*3-2:xi_id*3,:)), normals_storage.x0(s,:)', [zeros(3,3*nN_m),linear_get_centroid(cont_face{1}, s)]);
    Drotback_alg(s,xi_id*3-2:xi_id*3,:) = linear_reverse_rotate_points(slave_storage{s}.R, linear_rotation_matrix(normals_storage.n0(s,:)', [zeros(3,3*nN_m),linear_get_n0(cont_face{1}, normals_storage, s)]),...
      normals_storage.x0(s,:)', [zeros(3,3*nN_m),linear_get_centroid(cont_face{1}, s)], slave_storage{s}.rot(:,xi_id), squeeze(Drot_alg(s,xi_id*3-2:xi_id*3,:)));
    Dclip_alg(s,xi_id*3-2:xi_id*3,:) = linear_reverse_rotate_points(slave_storage{s}.R, linear_rotation_matrix(normals_storage.n0(s,:)', [zeros(3,3*nN_m),linear_get_n0(cont_face{1}, normals_storage, s)]),...
      normals_storage.x0(s,:)', [zeros(3,3*nN_m),linear_get_centroid(cont_face{1}, s)], clips_storage{s,1}.rot_vert(:,xi_id), squeeze(Drot_alg(s,xi_id*3-2:xi_id*3,:)));
  end
  for i=1:nN_m+nN_s
      if i <= nN_m
        for d=1:3
          cont_face_2 = cont_face{2};
          cont_face_2.coo(i,d) = cont_face_2.coo(i,d) + epsilon;
          [mort_D_2, mort_M_2, weight_gap_2, clips_storage_2, slave_storage_2, master_storage] = get_contact_dual_mortar(cont_face{1}, cont_face_2, normals_storage);
          Dproj_lim(s,:,d + (i-1)*3) = (reshape(slave_storage_2{s}.proj,[],1)-reshape(slave_storage{s}.proj,[],1))/epsilon;
          Drot_lim(s,:,d + (i-1)*3) = (reshape(slave_storage_2{s}.rot,[],1)-reshape(slave_storage{s}.rot,[],1))/epsilon;
          clip2 = reshape((slave_storage_2{s}.R' * (slave_storage_2{s}.rot-normals_storage.x0')) + normals_storage.x0',[],1);
          clip1 = reshape((slave_storage{s}.R' * (slave_storage{s}.rot-normals_storage.x0')) + normals_storage.x0',[],1);
          Drotback_lim(s,:,d + (i-1)*3) = (clip2-clip1)/epsilon;
          Dclip_lim(s,:,d + (i-1)*3) = (reshape(clips_storage_2{s,1}.vert,[],1)-reshape(clips_storage{s,1}.vert,[],1))/epsilon;
        end
      else
        for d=1:3
          cont_face_2 = cont_face{1};
          cont_face_2.coo(i-nN_m,d) = cont_face_2.coo(i-nN_m,d) + epsilon;
          normals_storage_2 = get_normals_centroids(cont_face_2);
          [mort_D_2, mort_M_2, weight_gap_2, clips_storage_2, slave_storage_2, master_storage] = get_contact_dual_mortar(cont_face_2, cont_face{2}, normals_storage_2);
          Dproj_lim(s,:,d + (i-1)*3) = (reshape(slave_storage_2{s}.proj,[],1)-reshape(slave_storage{s}.proj,[],1))/epsilon;
          Drot_lim(s,:,d + (i-1)*3) = (reshape(slave_storage_2{s}.rot,[],1)-reshape(slave_storage{s}.rot,[],1))/epsilon;
          clip2 = reshape((slave_storage_2{s}.R' * (slave_storage_2{s}.rot-normals_storage_2.x0')) + normals_storage_2.x0',[],1);
          clip1 = reshape((slave_storage{s}.R' * (slave_storage{s}.rot-normals_storage.x0')) + normals_storage.x0',[],1);
          Drotback_lim(s,:,d + (i-1)*3) = (clip2-clip1)/epsilon;
          Dclip_lim(s,:,d + (i-1)*3) = (reshape(clips_storage_2{s,1}.vert,[],1)-reshape(clips_storage{s,1}.vert,[],1))/epsilon;
          %[reshape(clips_storage_2{s,1}.vert,[],1)-clip2, reshape(clips_storage{s,1}.vert,[],1)-clip1]
        end
      end
   end
end
fprintf('Slave projections: %g \n', max(reshape(abs(Dproj_lim-Dproj_alg),[],1)));
fprintf('Slave rotated projections: %g \n', max(reshape(abs(Drot_lim-Drot_alg),[],1)));
fprintf('Slave rotated back projections: %g \n', max(reshape(abs(Drotback_lim-Drotback_alg),[],1)));
fprintf('Slave rotated back clip: %g \n', max(reshape(abs(Dclip_lim-Dclip_alg),[],1)));
%}
%{
% ----------------- clipping ----------------------------------- %
[mort_D, mort_M, weight_gap, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cont_face{1}, cont_face{2}, normals_storage);
ncells = size(clips_storage{1,1}.vert,2);
Dvert_alg = zeros(cont_face{1}.info.nele,cont_face{2}.info.nele,3*ncells,3*(nN_m+nN_s));
Dvert_lim = zeros(size(Dvert_alg));
for s=1:cont_face{1}.info.nele
  for m=1:cont_face{2}.info.nele
    V = linear_rot_polyg_vert(cont_face{1}.coo(cont_face{1}.nod(s,:),:)', cont_face{2}.coo(cont_face{2}.nod(m,:),:)', ...
      cont_face{1}.nod(s,:), cont_face{2}.nod(m,:), normals_storage.x0(s,:)', normals_storage.n0(s,:)',...
      [zeros(3,3*nN_m),linear_get_n0(cont_face{1}, normals_storage, s)],[zeros(3,3*nN_m), linear_get_centroid(cont_face{1}, s)], ...
      clips_storage{s,m}.orig, slave_storage{s}.proj, master_storage{s,m}.proj,...
      slave_storage{s}.rot, master_storage{s,m}.rot, nN_m, slave_storage{s}.R,...
      linear_rotation_matrix(normals_storage.n0(s,:)', [zeros(3,3*nN_m),linear_get_n0(cont_face{1}, normals_storage, s)]), clips_storage{s,m}.rot_vert);
    for cell=1:ncells
      Dvert_alg(s,m,cell*3-2:cell*3,:) = V{cell,2};
    end
    for i=1:nN_m+nN_s
      if i <= nN_m
        for d=1:3
          cont_face_2 = cont_face{2};
          cont_face_2.coo(i,d) = cont_face_2.coo(i,d) + epsilon;
          [mort_D_2, mort_M_2, weight_gap_2, clips_storage_2, slave_storage, master_storage] = get_contact_dual_mortar(cont_face{1}, cont_face_2, normals_storage);
          Dvert_lim(s,m,:,d + (i-1)*3) = (reshape(clips_storage_2{s,m}.vert,[],1)-reshape(clips_storage{s,m}.vert,[],1))/epsilon;
        end
      else
        for d=1:3
          cont_face_2 = cont_face{1};
          cont_face_2.coo(i-nN_m,d) = cont_face_2.coo(i-nN_m,d) + epsilon;
          normals_storage_2 = get_normals_centroids(cont_face_2);
          [mort_D_2, mort_M_2, weight_gap_2, clips_storage_2, slave_storage, master_storage] = get_contact_dual_mortar(cont_face_2, cont_face{2}, normals_storage_2);
          Dvert_lim(s,m,:,d + (i-1)*3) = (reshape(clips_storage_2{s,m}.vert,[],1)-reshape(clips_storage{s,m}.vert,[],1))/epsilon;
        end
      end
    end
  end
end
fprintf('Polygon vertices: %g \n', max(reshape(abs(Dvert_lim-Dvert_alg),[],1)));
%}


% ----- mortar matrices and weighted gap ---------------
[mort_D, mort_M, weight_gap, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cont_face{1}, cont_face{2}, normals_storage);
[mort_Dalg, mort_Malg, weight_gap_alg] = linear_dual_mortar_fem(cont_face{1}, cont_face{2}, clips_storage, normals_storage, slave_storage, master_storage);

z = ones(3*nN_s,1);

DDz_alg = squeeze(pagemtimes(mort_Dalg,z));
DMTz_alg = squeeze(pagemtimes(permute(mort_Malg, [2, 1, 3]),z));
DDz_lim = zeros(size(DDz_alg)); DMTz_lim = zeros(size(DMTz_alg));
Dweight_gap_lim = zeros(size(weight_gap_alg));
Dweight_gap_alg = weight_gap_alg;
for i=1:nN_m+nN_s
  if i <= nN_m
    for d=1:3
      cont_face_2 = cont_face{2};
      cont_face_2.coo(i,d) = cont_face_2.coo(i,d) + epsilon;
      [mort_D_2, mort_M_2, weight_gap_2, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cont_face{1}, cont_face_2, normals_storage);
      DDz_lim(:,d+(i-1)*3) = (mort_D_2*z-mort_D*z)/epsilon;
      DMTz_lim(:,d+(i-1)*3) = (mort_M_2'*z-mort_M'*z)/epsilon;
      Dweight_gap_lim(:,d+(i-1)*3) = (weight_gap_2-weight_gap)/epsilon;
    end
  else
    for d=1:3
      cont_face_2 = cont_face{1};
      cont_face_2.coo(i-nN_m,d) = cont_face_2.coo(i-nN_m,d) + epsilon;
      normals_storage_2 = get_normals_centroids(cont_face_2);
      [mort_D_2, mort_M_2, weight_gap_2, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cont_face_2, cont_face{2}, normals_storage_2);
      DDz_lim(:,d+(i-1)*3) = (mort_D_2*z-mort_D*z)/epsilon;
      DMTz_lim(:,d+(i-1)*3) = (mort_M_2'*z-mort_M'*z)/epsilon;
      Dweight_gap_lim(:,d+(i-1)*3) = (weight_gap_2-weight_gap)/epsilon;
    end
  end
end
fprintf('Difference between dir derivative calculated algebraically and with finite differences for:\n');
fprintf('Mortar matrix D times z: %g\n', max(reshape(abs(DDz_lim-DDz_alg),[],1)));
fprintf('Transpose of mortar matrix M times z: %g\n', max(reshape(abs(DMTz_lim-DMTz_alg),[],1)));
fprintf('weighted gap: %g\n', max(reshape(abs(Dweight_gap_lim-Dweight_gap_alg),[],1)));

%{
figure;
subplot(3,1,1);
plot(abs(DDz_lim(:)-DDz_alg(:)), 'b');
title('Rozdíl: \mathbf{D} \cdot \vec{z} (FD vs algebraicky)');
ylabel('Chyba');
grid on;

subplot(3,1,2);
plot(abs(DMTz_lim(:)-DMTz_alg(:)), 'r');
title('Rozdíl: \mathbf{M}^\top \cdot \vec{z} (FD vs algebraicky)');
ylabel('Chyba');
grid on;
%}

max_abs_diff_DDz = max(abs(DDz_lim - DDz_alg), [], 1);
max_abs_diff_DMTz = max(abs(DMTz_lim - DMTz_alg), [], 1);
max_abs_diff_wgap = max(abs(Dweight_gap_lim - Dweight_gap_alg), [], 1);
max_abs_diff_T1 = max(abs(Dt1_alg-Dt1_lim), [], 1);
max_abs_diff_T2 = max(abs(Dt2_alg-Dt2_lim), [], 1);


figure(1);           % vytvoří nebo přepne na figure s ID 1
set(gcf, 'Color', 'w');  % nastaví bílou barvu pozadí aktuálního okna
% Horní subplot – DDz rozdíl
subplot(2,1,1);
semilogy(max_abs_diff_DDz, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', [0.2 0.4 0.8], 'MarkerFaceColor', [0.2 0.4 0.8]);
grid on;
set(gca, 'FontSize', 11);  % skryje popisky X osy

% Dolní subplot – DMTz rozdíl
subplot(2,1,2);
semilogy(max_abs_diff_DMTz, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', [0.8 0.3 0.2], 'MarkerFaceColor', [0.8 0.3 0.2]);
grid on;
set(gca, 'FontSize', 11);

figure(2);           % vytvoří nebo přepne na figure s ID 1
set(gcf, 'Color', 'w');  % nastaví bílou barvu pozadí aktuálního okna
semilogy(max_abs_diff_wgap, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', [0.2 0.4 0.8], 'MarkerFaceColor', [0.2 0.4 0.8]);
grid on;
set(gca, 'FontSize', 11);

figure(3);           % vytvoří nebo přepne na figure s ID 1
set(gcf, 'Color', 'w');  % nastaví bílou barvu pozadí aktuálního okna
% Horní subplot – DDz rozdíl
subplot(2,1,1);
semilogy(max_abs_diff_T1, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', [0.2 0.4 0.8], 'MarkerFaceColor', [0.2 0.4 0.8]);
grid on;
set(gca, 'FontSize', 11);  % skryje popisky X osy

% Dolní subplot – DMTz rozdíl
subplot(2,1,2);
semilogy(max_abs_diff_T2, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', [0.8 0.3 0.2], 'MarkerFaceColor', [0.8 0.3 0.2]);
grid on;
set(gca, 'FontSize', 11);

