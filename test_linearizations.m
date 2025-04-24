clear
%addpath(genpath('..'));
addpath("core");
addpath("plot");
addpath("core_matfem");
addpath(genpath("linearizations"));
addpath("clipper2");
addpath("mtimesx");


mesh = cell(2,1);
etype = cell(2,1);
etype{1} = 'hexa8';
etype{2} = 'hexa8';

% Create reference configurations X^(1,2), t=0
mesh{1} = mesh_generator([6,6,1], [3,3,2],etype{1}); 
mesh{2} = mesh_generator([9,9,1], [3,3,2],etype{2});

mesh{1}.nod.coo=mesh{1}.nod.coo + [0,0,3];
mesh{1}.nod.coo=mesh{1}.nod.coo + [0,2,0];
mesh{1}.nod.coo=mesh{1}.nod.coo + [2,0,0];

% Changes to the geometry
mesh{1}.nod.coo(:, 3) = mesh{1}.nod.coo(:, 3) - mesh{1}.nod.coo(:, 1) / 5;
%mesh{1}.nod.coo(mesh{1}.bou{5}.nod(1,1), :) = mesh{1}.nod.coo(mesh{1}.bou{5}.nod(1,1), :) - 0.5;
%mesh{1}.nod.coo(mesh{1}.bou{5}.nod(1,2), :) = mesh{1}.nod.coo(mesh{1}.bou{5}.nod(1,2), :) - 0.5;
%mesh{2}.nod.coo(mesh{2}.bou{6}.nod(1,2), :) = mesh{2}.nod.coo(mesh{2}.bou{6}.nod(1,2), :) - 0.5;

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
for i=8:14
  for j=1:4
    cont_face_2 = cont_face;
    delt_d_c = zeros(3*nN_m +3*nN_s, 1);
    delt_di = 1e-6*rand(3,1);
    random_id = randi(nN_m +nN_s);
    delt_d_c(random_id*3-2:random_id*3) = delt_di;
    delt_d_s = reshape(delt_d_c(3*nN_m+1:end),3,nN_s)';
    delt_d_m = reshape(delt_d_c(1:3*nN_m),3,nN_m)';
  
    
    delt_Dz_mat = squeeze(pagemtimes(mort_Dalg,z))*delt_d_c;
    delt_Mz_mat = squeeze(pagemtimes(permute(mort_Malg, [2, 1, 3]),z))*delt_d_c;
    delt_wgap_mat = weight_gap_alg*delt_d_c;
    epsilon = 10^-i;
    cont_face_2{1}.coo = cont_face{1}.coo + epsilon*delt_d_s;   % f(x + epsilon*delt_d)
    cont_face_2{2}.coo = cont_face{2}.coo + epsilon*delt_d_m;
    
    normals_storage_2 = get_normals_centroids(cont_face_2{1});
    [mort_D_2, mort_M_2, weight_gap_2, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cont_face_2{1}, cont_face_2{2}, normals_storage_2);
    delt_Dz = (mort_D_2*z-mort_D*z)/epsilon;
    delt_Mz = (mort_M_2'*z-mort_M'*z)/epsilon;
    delt_wgap = (weight_gap_2-weight_gap)/epsilon;
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
  end

end

%{
eps_vals = 10.^(-2:-1:-12);
errors_D = zeros(size(eps_vals));
errors_M = zeros(size(eps_vals));

for i = 1:length(eps_vals)
    epsilon = eps_vals(i);
    cont_face_2{1}.coo = cont_face{1}.coo + epsilon*delt_d_s;   % f(x + epsilon*delt_d)
    cont_face_2{2}.coo = cont_face{2}.coo + epsilon*delt_d_m;
  
    normals_storage_2 = get_normals_centroids(cont_face_2{1});
    [mort_D_2, mort_M_2, weight_gap, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cont_face_2{1}, cont_face_2{2}, normals_storage_2);
    delt_Dz = (mort_D_2*z-mort_D*z)/epsilon;
    delt_Mz = (mort_M_2'*z-mort_M'*z)/epsilon;
    errors_D(i) = norm(delt_Dz - delt_Dz_mat); % or relative error
    errors_M(i) = norm(delt_Mz - delt_Mz_mat); % or relative error
end

loglog(eps_vals, errors_D, '-o');
xlabel('Epsilon');
ylabel('Error');
title('Convergence of Finite Difference Approximation for matrix D');
grid on;
%}

%{
loglog(eps_vals, errors_M, '-o');
xlabel('Epsilon');
ylabel('Error');
title('Convergence of Finite Difference Approximation for matrix M');
grid on;
%}

