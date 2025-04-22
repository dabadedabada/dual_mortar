clear
%addpath(genpath('..'));
addpath("core");
addpath("plot");
addpath("core_matfem");
addpath(genpath("linearizations"));
addpath("clipper2");

mesh = cell(2,1);
etype = cell(2,1);
etype{1} = 'hexa8';
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

% PLOT ------------
%plot_meshes(mesh);

disc = cell(2,1);
bypr = cell(2,1);

FE = cell(2,1);
FE{1} = fe_init(etype{1},2);
FE{2} = fe_init(etype{2},2);

consts = cell(2,1);
consts{1} = struct('Young', 1e3, 'Poisson', 0.3, 'density', 1e1, 'lambda', [], 'mu', []);
consts{2} = struct('Young', 2e3, 'Poisson', 0.4, 'density', 2e1, 'lambda', [], 'mu', []);

for i=1:2
  consts{i}.lambda = consts{i}.Young*consts{i}.Poisson/((1+consts{i}.Poisson)*(1-2*consts{i}.Poisson));
  consts{i}.mu = consts{i}.Young/(2*(1+consts{i}.Poisson));
end

times = 0:0.1:1;

problem = cell(2,1);
problem{1} = struct('material_model', 'kirchhoff', 'volume_force', [0,0,0], 'variables', [], 'boundary', struct('dirichlet', []));
problem{1}.boundary.dirichlet = {...
  6,'0*[X,Y,Z]+[ 1,0,0]*T','add'};
problem{2} = struct('material_model', 'kirchhoff', 'volume_force', [0,0,0], 'variables', [], 'boundary', struct('dirichlet', []));
problem{2}.boundary.dirichlet = {...
  1,'0*[X,Y,Z]+[ 0,0,0]','add'; ...
  2,'0*[X,Y,Z]+[ 0,0,0]','add'; ...
  3,'0*[X,Y,Z]+[ 0,0,0]','add'; ...
  4,'0*[X,Y,Z]+[ 0,0,0]','add'};

u = cell(2,1);
u{1} = zeros(3*mesh{1}.info.nodcount,1);
u{2} = zeros(3*mesh{2}.info.nodcount,1);

Bc = cell(2,1);

for time_ind=1:size(times,2)-1
  time=times(time_ind);

  [disc{1}, bypr{1}] = assemble_hyperelasticity(FE{1}, mesh{1}, consts{1}, problem{1}, u{1});
  [disc{2}, bypr{2}] = assemble_hyperelasticity(FE{2}, mesh{2}, consts{2}, problem{2}, u{2});

  %Bc{1} = set_boundary_conditions(mesh{1}, problem{1}, time, time_ind, u{1} );
  %Bc{2} = set_boundary_conditions(mesh{2}, problem{2}, time, time_ind, u{2} );

  % Init contact faces
  cont_face = cell(2,1);
  cont_face{1} = init_cont_face(mesh{1}, 5);
  cont_face{2} = init_cont_face(mesh{2}, 6);

  % semi-smooth Newton
  % calculate linearizations for B_co(d)=[0;-M^T;D^T] given displacements d
  % Popp diss (4.23)

  % Normals, centroids and tangents
  normals_storage = get_normals_centroids(cont_face{1});

  % Algorithm 1
  [mort_D, mort_M, weight_gap, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cont_face{1}, cont_face{2}, normals_storage);

  % Algorithm 1 linearization
  %Ctilde = linear_dual_mortar_fem(cont_face{1}, cont_face{2}, clips_storage, normals_storage);

end





