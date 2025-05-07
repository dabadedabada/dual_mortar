clear
%addpath(genpath('..'));
addpath("core");
addpath("plot");
addpath("core_matfem");
addpath(genpath("linearizations"));
addpath("clipper2");

mesh = cell(2,1);
etype = cell(2,1);
etype{1} = 'tetr4';
etype{2} = 'hexa8';

% Create reference configurations X^(1,2), t=0
mesh{1} = generate_cylinder_mesh(1,9, [10,10,10],etype{1});
%mesh{1} = mesh_generator([3,4,2], [5,5,5],etype{1}); 
mesh{2} = mesh_generator([9,9,1], [4,3,2],etype{2});

mesh{1}.nod.coo=mesh{1}.nod.coo + [0,0,3];
%mesh{1}.nod.coo=mesh{1}.nod.coo + [0,2,0];
%mesh{1}.nod.coo=mesh{1}.nod.coo + [2,0,0];

% Changes to the geometry
%mesh{1}.nod.coo(:, 3) = mesh{1}.nod.coo(:, 3) - mesh{1}.nod.coo(:, 1) / 5;
%mesh{1}.nod.coo(mesh{1}.bou{5}.nod(1,1), :) = mesh{1}.nod.coo(mesh{1}.bou{5}.nod(1,1), :) - 0.5;
%mesh{1}.nod.coo(mesh{1}.bou{5}.nod(1,2), :) = mesh{1}.nod.coo(mesh{1}.bou{5}.nod(1,2), :) - 0.5;
%mesh{2}.nod.coo(mesh{2}.bou{6}.nod(1,2), :) = mesh{2}.nod.coo(mesh{2}.bou{6}.nod(1,2), :) - 0.5;

plot_meshes(mesh);