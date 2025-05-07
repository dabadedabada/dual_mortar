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
mesh{1} = mesh_generator([3,3,2], [1,1,1],etype{1});
mesh{2} = mesh_generator([12,12,1], [1,1,1],etype{2});

%mesh{1}.nod.coo=mesh{1}.nod.coo + [0,0,1.1];
%mesh{1}.nod.coo=mesh{1}.nod.coo + [0.5,1.5,0];

mesh{1}.nod.coo = mesh{1}.nod.coo + [-1.5,-1.5,0 ];
mesh{2}.nod.coo = mesh{2}.nod.coo + [-6.0,-6.0,0 ];
if true
  mesh{1}.nod.coo = mesh{1}.nod.coo + [0,0,1.01];
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

% -------- test linearizations of normals, tangents, n0 and x0 ------------
%test_linear_normals_centroids(cont_face{1}, cont_face{2}, delt_d_s, delt_d_m); 

% Algorithm 1
[mort_D, mort_M, weight_gap, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cont_face{1}, cont_face{2}, normals_storage);
[mort_Dalg, mort_Malg, weight_gap_alg] = linear_dual_mortar_fem(cont_face{1}, cont_face{2}, clips_storage, normals_storage, slave_storage, master_storage);



z = ones(3*nN_s,1);
err = 10^-4;
epsilon = 10^-8;
for i=1:3*(nN_s+nN_m)
    cont_face_2 = cont_face;
    delt_d_c = zeros(3*nN_m +3*nN_s, 1);
    %random_id = randi(nN_m +nN_s);
    delt_d_c(i) = 0.001;
    delt_d_s = reshape(delt_d_c(3*nN_m+1:end),3,nN_s)';
    delt_d_m = reshape(delt_d_c(1:3*nN_m),3,nN_m)';
  
    delt_Dz_mat = squeeze(pagemtimes(mort_Dalg,z))*delt_d_c;
    delt_Mz_mat = squeeze(pagemtimes(permute(mort_Malg, [2, 1, 3]),z))*delt_d_c;
    delt_wgap_mat = weight_gap_alg*delt_d_c;

    cont_face_2{1}.coo = cont_face{1}.coo + epsilon*delt_d_s;   % f(x + epsilon*delt_d)
    cont_face_2{2}.coo = cont_face{2}.coo + epsilon*delt_d_m;
    
    normals_storage_2 = get_normals_centroids(cont_face_2{1});

    fprintf('DOF: %d\n', i);
    
    if (i > 3*nN_m)
      for s=1:cont_face{1}.info.nele
        for xi_id=1:cont_face{1}.info.nN_ele
          Nhatje = linear_get_nodal_element_normals(cont_face{1}, xi_id, s);
          delt_nhatje_mat = Nhatje*reshape(delt_d_s',[3*nN_s, 1]);
          delt_nhatje = (normals_storage_2.ele_normals{s}(xi_id,:)-normals_storage.ele_normals{s}(xi_id,:))/epsilon;
          if (norm(delt_nhatje-delt_nhatje_mat) < err)
             fprintf('For epsilon = %g: nhatje  s: %d, xi: %d         - Success - with difference of dir derivatives %g \n', epsilon, s, xi_id, norm(delt_nhatje-delt_nhatje_mat));
          else
             fprintf('For epsilon = %g: nhatje  s: %d, xi: %d         - Fail    - with difference of dir derivatives %g \n', epsilon, s, xi_id, norm(delt_nhatje-delt_nhatje_mat));
          end

        end
        %{
        N0 = linear_get_n0(cont_face{1}, normals_storage, s);
        delt_n0_mat = N0*reshape(delt_d_s',[3*nN_s, 1]);
        delt_n0 = (normals_storage_2.n0(s,:)-normals_storage.n0(s,:))/epsilon;
        if (norm(delt_n0-delt_n0_mat) < err)
           fprintf('For epsilon = %g: n0  s: %d         - Success - with difference of dir derivatives %g \n', epsilon, s, norm(delt_n0-delt_n0_mat));
        else
           fprintf('For epsilon = %g: n0  s: %d         - Fail    - with difference of dir derivatives %g \n', epsilon, s, norm(delt_n0-delt_n0_mat));
        end
        %}

      end
      
      
    end
    

    %{
    [mort_D_2, mort_M_2, weight_gap_2, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cont_face_2{1}, cont_face_2{2}, normals_storage_2);
    delt_Dz = (mort_D_2*z-mort_D*z)/epsilon;
    delt_Mz = (mort_M_2'*z-mort_M'*z)/epsilon;
    delt_wgap = (weight_gap_2-weight_gap)/epsilon;
    fprintf('DOF: %d\n', i);
    if (norm(delt_Dz-delt_Dz_mat) < err)
       fprintf('For epsilon = %g: Dz           - Success - with difference of dir derivatives %g \n', epsilon, norm(delt_Dz-delt_Dz_mat));
    else
       fprintf('For epsilon = %g: Dz           - Fail    - with difference of dir derivatives %g \n', epsilon, norm(delt_Dz-delt_Dz_mat));
    end
    if (norm(delt_Mz-delt_Mz_mat) < err)
       fprintf('For epsilon = %g: Mz           - Success - with difference of dir derivatives %g \n', epsilon, norm(delt_Mz-delt_Mz_mat));
    else
       fprintf('For epsilon = %g: Mz           - Fail    - with difference of dir derivatives %g \n', epsilon, norm(delt_Mz-delt_Mz_mat));
    end
    if (norm(delt_wgap-delt_wgap_mat) < err)
       fprintf('For epsilon = %g: weighted_gap - Success - with difference of dir derivatives %g \n', epsilon, norm(delt_wgap-delt_wgap_mat));
    else
       fprintf('For epsilon = %g: weighted_gap - Fail    - with difference of dir derivatives %g \n', epsilon, norm(delt_wgap-delt_wgap_mat));
    end
    %}
end


