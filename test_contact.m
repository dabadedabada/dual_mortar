clear
%addpath(genpath('..'));
addpath("core");
%addpath("plot");
addpath("core_matfem");
addpath("linearizations");
addpath("clipper2");

mesh = cell(2,1);
etype = cell(2,1);
etype{1} = 'hexa8';
etype{2} = 'tetr4';

% Create reference configurations X^(1,2), t=0
mesh{1} = mesh_generator([6,6,1], [3,3,2],etype{1}); 
mesh{2} = mesh_generator([9,9,1], [3,3,2],etype{2});

mesh{1}.nod.coo=mesh{1}.nod.coo + [0,0,0.9];
mesh{1}.nod.coo=mesh{1}.nod.coo + [0.5,1.5,0];

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
consts{1} = struct('Young', 1e3, 'Poisson', 0.3, 'density', 1e-1, 'lambda', [], 'mu', []);
consts{2} = struct('Young', 2e3, 'Poisson', 0.4, 'density', 2e-1, 'lambda', [], 'mu', []);

for i=1:2
  consts{i}.lambda = consts{i}.Young*consts{i}.Poisson/((1+consts{i}.Poisson)*(1-2*consts{i}.Poisson));
  consts{i}.mu = consts{i}.Young/(2*(1+consts{i}.Poisson));
end

problem = cell(2,1);
problem{1} = struct('material_model', 'kirchhoff', 'volume_force', [0,0,0], 'variables', [], 'boundary', struct('dirichlet', []));
problem{1}.boundary.dirichlet = {...
  6,'0*[X,Y,Z]+[ 1.5,0,0]*T+[0,0,-0.1]','enforce'};
problem{2} = struct('material_model', 'kirchhoff', 'volume_force', [0,0,0], 'variables', [], 'boundary', struct('dirichlet', []));
problem{2}.boundary.dirichlet = {...
  1,'0*[X,Y,Z]+[ 0,0,0]','enforce'; ...
  2,'0*[X,Y,Z]+[ 0,0,0]','enforce'; ...
  3,'0*[X,Y,Z]+[ 0,0,0]','enforce'; ...
  4,'0*[X,Y,Z]+[ 0,0,0]','enforce'};

u = cell(2,1);
u{1} = zeros(3*mesh{1}.info.nodcount,1);
u{2} = zeros(3*mesh{2}.info.nodcount,1);

BC = cell(2,1);

timestepping.steps = 10;
timestepping.times = 0:1/timestepping.steps:1;

ensight_name = 'slide';
exports = {'Displacement','nod.vec';'Stress','ele.tensym'};
export_ensight_geometry(    ensight_name, mesh{1}, timestepping, exports, [min(mesh{1}.info.extent(:,1),mesh{2}.info.extent(:,1)), max(mesh{1}.info.extent(:,2),mesh{2}.info.extent(:,2))]);
export_ensight_geometry_add(ensight_name, mesh{2}, 1);
export_ensight_result(      ensight_name, reshape(zeros(mesh{1}.info.nodcount,3)',[],1), struct(  'type','nod.vec'   ,'eletype',mesh{1}.ele.type,'name','Displacement'),'000'    );
export_ensight_result_add(  ensight_name, reshape(zeros(mesh{2}.info.nodcount,3)',[],1), struct(  'type','nod.vec'   ,'eletype',mesh{2}.ele.type,'name','Displacement'),'000',3,1);
export_ensight_result(      ensight_name, reshape(zeros(mesh{1}.info.elecount,6)',[],1), struct(  'type','ele.tensym','eletype',mesh{1}.ele.type,'name','Stress'      ),'000'    );
export_ensight_result_add(  ensight_name, reshape(zeros(mesh{2}.info.elecount,6)',[],1), struct(  'type','ele.tensym','eletype',mesh{2}.ele.type,'name','Stress'      ),'000',3,1);

for n = 1:timestepping.steps
  t  =   timestepping.times(n+1);
  dt = t-timestepping.times(n); 

  [disc{1}, bypr{1}] = assemble_hyperelasticity(FE{1}, mesh{1}, consts{1}, problem{1}, u{1});
  [disc{2}, bypr{2}] = assemble_hyperelasticity(FE{2}, mesh{2}, consts{2}, problem{2}, u{2});

  BC{1} = set_boundary_conditions(mesh{1}, problem{1}, t, n, u{1} );
  BC{2} = set_boundary_conditions(mesh{2}, problem{2}, t, n, u{2} );

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

  DK = blkdiag( disc{1}.Kc+disc{1}.Ks+(dt^-2)*disc{1}.M,  disc{2}.Kc+disc{2}.Ks+(dt^-2)*disc{2}.M);
  r  = [disc{1}.fext+disc{1}.fint+(dt^-2)*disc{1}.M*u{1}; disc{2}.fext+disc{2}.fint+(dt^-2)*disc{2}.M*u{2}];
  
  options = optimoptions('quadprog','Display','none');

  tmp  = quadprog(...
    DK, -r,...
    [],[],...
    blkdiag(BC{1}.Bd{1},BC{2}.Bd{1}), [BC{1}.Bd{1}*BC{1}.U0;BC{2}.Bd{1}*BC{2}.U0],...
    [],[],[],options);
  u{1} = tmp(             1:size(u{1},1));
  u{2} = tmp(size(u{1},1)+1:size(u{1},1)+size(u{2},1));

  export_ensight_result(     ensight_name, reshape(u{1}              ,[],1), struct(  'type','nod.vec'   ,'eletype',mesh{1}.ele.type,'name','Displacement'),sprintf('%03d',n)    );
  export_ensight_result_add( ensight_name, reshape(u{2}              ,[],1), struct(  'type','nod.vec'   ,'eletype',mesh{2}.ele.type,'name','Displacement'),sprintf('%03d',n),3,1);
  export_ensight_result(     ensight_name, reshape(bypr{1}.ele.sigma',[],1), struct(  'type','ele.tensym','eletype',mesh{1}.ele.type,'name','Stress'      ),sprintf('%03d',n)    );
  export_ensight_result_add( ensight_name, reshape(bypr{2}.ele.sigma',[],1), struct(  'type','ele.tensym','eletype',mesh{2}.ele.type,'name','Stress'      ),sprintf('%03d',n),3,1);


end





