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
etype{1} = 'tetr4';
etype{2} = 'tetr4';

% Create reference configurations X^(1,2), t=0
mesh{1} = mesh_generator([6,6,1], [3,3,2],etype{1}); 
mesh{2} = mesh_generator([9,9,1], [3,3,2],etype{2});

mesh{1}.nod.coo=mesh{1}.nod.coo + [0,0,3];
mesh{1}.nod.coo=mesh{1}.nod.coo + [0,2,0];
mesh{1}.nod.coo=mesh{1}.nod.coo + [2,0,0];

% Changes to the geometry
%mesh{1}.nod.coo(:, 3) = mesh{1}.nod.coo(:, 3) - mesh{1}.nod.coo(:, 1) / 5;
%mesh{1}.nod.coo(mesh{1}.bou{5}.nod(1,1), :) = mesh{1}.nod.coo(mesh{1}.bou{5}.nod(1,1), :) - 0.5;
%mesh{1}.nod.coo(mesh{1}.bou{5}.nod(1,2), :) = mesh{1}.nod.coo(mesh{1}.bou{5}.nod(1,2), :) - 0.5;
%mesh{2}.nod.coo(mesh{2}.bou{6}.nod(1,2), :) = mesh{2}.nod.coo(mesh{2}.bou{6}.nod(1,2), :) - 0.5;

% Init contact faces
cont_face = cell(2,1);
cont_face{1} = init_cont_face(mesh{1}, 5);
cont_face{2} = init_cont_face(mesh{2}, 6);

% Normals, centroids and tangents
normals_storage = get_normals_centroids(cont_face{1});

nN_s = cont_face{1}.info.nN;
d_s = rand(nN_s,3);
delt_d_s = rand(nN_s,3);

% -------- test linearizations of normals, tangents, n0 and x0 ------------
%test_linear_normals_centroids(cont_face{1}, d_s, delt_d_s); 

% Algorithm 1
[mort_D, mort_M, weight_gap, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cont_face{1}, cont_face{2}, normals_storage);

Ctilde = linear_dual_mortar_fem(cont_face{1}, cont_face{2}, clips_storage, normals_storage, slave_storage, master_storage);


