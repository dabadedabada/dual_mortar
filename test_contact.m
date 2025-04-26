clear
%addpath(genpath('..'));
%addpath("plot");
addpath('core_matfem');
addpath('core');
addpath('clipper2');
%addpath('../../../../programming/matlab/matfem_marta/core');
%addpath("linearizations");

mesh  = cell(2,1);
etype = cell(2,1);
etype{1} = 'hexa8';
etype{2} = 'tetr4';

% Create reference configurations X^(1,2), t=0
mesh{1} = mesh_generator([5,5,2], [6,6,2],etype{1});
mesh{2} = mesh_generator([9,9,1], [8,8,2],etype{2});

mesh{1}.nod.coo=mesh{1}.nod.coo + [0,0,1.1];
mesh{1}.nod.coo=mesh{1}.nod.coo + [0.5,1.5,0];

% Changes to the geometry
%mesh{1}.nod.coo(:, 3) = mesh{1}.nod.coo(:, 3) - mesh{1}.nod.coo(:, 1) / 5;
%mesh{1}.nod.coo(mesh{1}.bou{5}.nod(1,1), :) = mesh{1}.nod.coo(mesh{1}.bou{5}.nod(1,1), :) - 0.5;
%mesh{1}.nod.coo(mesh{1}.bou{5}.nod(1,2), :) = mesh{1}.nod.coo(mesh{1}.bou{5}.nod(1,2), :) - 0.5;
%mesh{2}.nod.coo(mesh{2}.bou{6}.nod(1,2), :) = mesh{2}.nod.coo(mesh{2}.bou{6}.nod(1,2), :) - 0.5;

% PLOT ------------
%plot_meshes(mesh);

% Init contact faces
cont_face = cell(2,1);
cont_face{1} = init_cont_face(mesh{1}, 5);
cont_face{2} = init_cont_face(mesh{2}, 6);

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

compl_param = consts{1}.Young;

problem = cell(2,1);
problem{1} = struct('material_model', 'kirchhoff', 'volume_force', [0,0,-10], 'variables', [], 'boundary', struct('dirichlet', []));
problem{1}.boundary.dirichlet = {...
  6,'0*[X,Y,Z]+[ 1.5,0,0]*max(T-0.2,0)+[0,0,-2.0]*min(T,0.2)','enforce'};
problem{2} = struct('material_model', 'kirchhoff', 'volume_force', [0,0,-10], 'variables', [], 'boundary', struct('dirichlet', []));
problem{2}.boundary.dirichlet = {...
  1,'0*[X,Y,Z]+[ 0,0,0]','enforce'; ...
  2,'0*[X,Y,Z]+[ 0,0,0]','enforce'; ...
  3,'0*[X,Y,Z]+[ 0,0,0]','enforce'; ...
  4,'0*[X,Y,Z]+[ 0,0,0]','enforce'};
% problem{2}.boundary.dirichlet = {...
%   5,'0*[X,Y,Z]+[ 0,0,0]','enforce'};

u = cell(2,1); X = cell(2,1);
u{1} = zeros(3*mesh{1}.info.nodcount,1);  X{1} = reshape(mesh{1}.nod.coo',[],1);
u{2} = zeros(3*mesh{2}.info.nodcount,1);  X{2} = reshape(mesh{2}.nod.coo',[],1);

BC = cell(2,1);

timestepping.steps = 20;
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

  fprintf('Time (%3i) %f \n',n,t);
  
  for c = 1:2
    fprintf(' contact %d\n',c);
    BC{1} = set_boundary_conditions(mesh{1}, problem{1}, t, n, u{1} );
    BC{2} = set_boundary_conditions(mesh{2}, problem{2}, t, n, u{2} );

    % Normals, centroids and tangents
    normals_storage = get_normals_centroids(cont_face{1});

    % Algorithm 1 - for reference configuration X
    [mort_D, mort_M, weight_gap, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cont_face{1}, cont_face{2}, normals_storage);
    [mort_Dalg, mort_Malg, weight_gap_alg] = linear_dual_mortar_fem(cont_face{1}, cont_face{2}, clips_storage, normals_storage, slave_storage, master_storage);
    tmp_normals = mat2cell(normals_storage.averaged_normals,ones(size(normals_storage.averaged_normals,1),1),3);
    tmp_normals = blkdiag(tmp_normals{:});
    mort_M = mort_D\mort_M;
    mort_D = mort_D\mort_D;
    ND = tmp_normals*mort_D; tmpmap = reshape([3*cont_face{1}.map.nod-2 3*cont_face{1}.map.nod-1 3*cont_face{1}.map.nod]',[],1); [i_,j_,v_] = find(ND); ND = sparse(i_,tmpmap(j_),v_,max(i_),mesh{1}.info.noddof*mesh{1}.info.nodcount);
    NM = tmp_normals*mort_M; tmpmap = reshape([3*cont_face{2}.map.nod-2 3*cont_face{2}.map.nod-1 3*cont_face{2}.map.nod]',[],1); [i_,j_,v_] = find(NM); NM = sparse(i_,tmpmap(j_),v_,max(i_),mesh{2}.info.noddof*mesh{2}.info.nodcount);
    gap = -[ND, -NM]*[X{1}; X{2}];


    inact_loc = 1:size(cont_face{1}.map.nod, 1);       % init inactive 0 2D indices
    activ_loc = [];                                    % init active 0   2D indices

    inact_glob = cont_face{1}.map.nod(inact_loc);      % init inactive 0 3D indices
    activ_glob = cont_face{1}.map.nod(activ_loc);      % init active 0   3D indices

    m_id_loc = 1:size(cont_face{2}.map.nod, 1);        % init set{M} 2D indices
    m_id_glob = cont_face{2}.map.nod;                  % init set{M} 3D indices

    z = zeros(3*size(cont_face{1}.map.nod, 1),1);       % init z0

    k = 1;
    newton_continue = true;
    while newton_continue
      [disc{1}, bypr{1}] = assemble_hyperelasticity(FE{1}, mesh{1}, consts{1}, problem{1}, u{1});
      [disc{2}, bypr{2}] = assemble_hyperelasticity(FE{2}, mesh{2}, consts{2}, problem{2}, u{2});

      DK = blkdiag( disc{1}.Kc+disc{1}.Ks+(0*dt^-2)*disc{1}.M,                disc{2}.Kc+disc{2}.Ks    +(0*dt^-2)*disc{2}.M);
      r  =       [(-disc{1}.fext+disc{1}.fint   +(0*dt^-2)*disc{1}.M*u{1}); (-disc{2}.fext+disc{2}.fint+(0*dt^-2)*disc{2}.M*u{2})];

      % z^(k+1) = 0 on inactive set
      %z(inact_cont) = 0; 
      idx_inact_loc = reshape([3*inact_loc-2; 3*inact_loc-1; 3*inact_loc], 1, []);
      idx_activ_loc = reshape([3*activ_loc-2; 3*activ_loc-1; 3*activ_loc], 1, []);
      idx_m_loc = reshape([3*m_id_loc-2; 3*m_id_loc-1; 3*m_id_loc], 1, []);

      Ctilde = [squeeze(pagemtimes(permute(mort_Malg, [2, 1, 3]),z));squeeze(pagemtimes(mort_Dalg,z))];
      Ctilde_MM = Ctilde(idx_m_loc,idx_m_loc);
      Ctilde_AM = Ctilde(idx_activ_loc+ 3*cont_face{2}.info.nN,idx_m_loc);
      Ctilde_AI = Ctilde(idx_activ_loc+ 3*cont_face{2}.info.nN,idx_inact_loc+ 3*cont_face{2}.info.nN);
      Ctilde_AA = Ctilde(idx_activ_loc+ 3*cont_face{2}.info.nN,idx_activ_loc+ 3*cont_face{2}.info.nN);
      Ftilde_AA = []; Ftilde_AI =[]; T_AA = [];
      if ~isempty(activ_loc)
        Ftilde_xi = zeros(length(activ_loc),3*cont_face{1}.info.nN); Ftilde_eta = zeros(length(activ_loc),3*cont_face{1}.info.nN);
        T_xi = zeros(length(activ_loc),3*cont_face{1}.info.nN); T_eta = zeros(length(activ_loc),3*cont_face{1}.info.nN);
        for i=1:length(activ_loc)
          [T1, T2] = linear_get_tangents(cont_face{1}, normals_storage, activ_loc(i));
          zi = z(idx_activ_loc(3*i-2:3*i));
          Ftilde_xi(i,:) = zi'*T1; Ftilde_eta(i,:) = zi'*T2;
          T_xi(i,idx_activ_loc(3*i-2:3*i)) = normals_storage.t1(activ_loc(i));
          T_eta(i,idx_activ_loc(3*i-2:3*i)) = normals_storage.t2(activ_loc(i));
        end
        Ftilde_A = [Ftilde_xi;Ftilde_eta];     
        Ftilde_AA = Ftilde_A(:,idx_activ_loc);   
        Ftilde_AI = Ftilde_A(:,idx_inact_loc);
        T_A = [T_xi;T_eta];
        T_AA = T_A(:,idx_activ_loc);  % toto je v jeho matici T_A
      end
      Stilde_AA = weight_gap_alg(activ_loc, idx_activ_loc + 3*cont_face{2}.info.nN);
      Stilde_AI = weight_gap_alg(activ_loc, idx_inact_loc + 3*cont_face{2}.info.nN);
      Mtilde_A =  weight_gap_alg(activ_loc, idx_m_loc);



      if false
        options = optimoptions('quadprog','Display','none');
        tmp  = quadprog(...
          (DK+DK')/2, -r,...
          [], [],...
          [blkdiag(BC{1}.Bd{1},BC{2}.Bd{1});...
          [ND, -NM]],...
          [BC{1}.Bd{1}*(BC{1}.U0-u{1});BC{2}.Bd{1}*(BC{2}.U0-u{2});...
          gap-[ND, -NM]*[u{1};u{2}]],...
          [],[],[],options);
      else
        B_ = [blkdiag(BC{1}.Bd{1},BC{2}.Bd{1})                     ;   [ND, -NM]];
        c_ = [BC{1}.Bd{1}*(BC{1}.U0-u{1}); BC{2}.Bd{1}*(BC{2}.U0-u{2});   gap-[ND, -NM]*[u{1};u{2}]];
        tmp = [DK B_'; B_ zeros(size(B_,1))]\[-r; c_];
      end
      du{1} = tmp(             1:size(u{1},1));
      du{2} = tmp(size(u{1},1)+1:size(u{1},1)+size(u{2},1));
      if k == 1,      newton_err = norm([du{1}; du{2}]);
      else,           newton_err = norm([du{1}; du{2}])/norm([u{1}; u{2}]);
      end
      fprintf('  newton %3d: %e\n',k,newton_err);
      newton_continue = (k <= 1) || (newton_err > 1e-6);
      u{1} = u{1} + du{1};
      u{2} = u{2} + du{2};
      export_ensight_result(     ensight_name, reshape(u{1}              ,[],1), struct(  'type','nod.vec'   ,'eletype',mesh{1}.ele.type,'name','Displacement'),sprintf('%03d',n+k-1)    );
      export_ensight_result_add( ensight_name, reshape(u{2}              ,[],1), struct(  'type','nod.vec'   ,'eletype',mesh{2}.ele.type,'name','Displacement'),sprintf('%03d',n+k-1),3,1);
      export_ensight_result(     ensight_name, reshape(bypr{1}.ele.sigma',[],1), struct(  'type','ele.tensym','eletype',mesh{1}.ele.type,'name','Stress'      ),sprintf('%03d',n+k-1)    );
      export_ensight_result_add( ensight_name, reshape(bypr{2}.ele.sigma',[],1), struct(  'type','ele.tensym','eletype',mesh{2}.ele.type,'name','Stress'      ),sprintf('%03d',n+k-1),3,1);

      % add d^k+1 to cont_face coordinates
      cont_face{1}.coo = cont_face{1}.coo + u{1}(cont_face{1}.map.nod); 
      cont_face{2}.coo = cont_face{2}.coo + u{2}(cont_face{2}.map.nod);

      % Normals, centroids and tangents for d^k+1
      normals_storage = get_normals_centroids(cont_face{1});
  
      % Algorithm 1 - for current configuration with d^k+1
      [mort_D, mort_M, weight_gap, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cont_face{1}, cont_face{2}, normals_storage);
      [mort_Dalg, mort_Malg, weight_gap_alg] = linear_dual_mortar_fem(cont_face{1}, cont_face{2}, clips_storage, normals_storage, slave_storage, master_storage);
      tmp_normals = mat2cell(normals_storage.averaged_normals,ones(size(normals_storage.averaged_normals,1),1),3);
      tmp_normals = blkdiag(tmp_normals{:});
      mort_M = mort_D\mort_M;
      mort_D = mort_D\mort_D;
      ND = tmp_normals*mort_D; tmpmap = reshape([3*cont_face{1}.map.nod-2 3*cont_face{1}.map.nod-1 3*cont_face{1}.map.nod]',[],1); [i_,j_,v_] = find(ND); ND = sparse(i_,tmpmap(j_),v_,max(i_),mesh{1}.info.noddof*mesh{1}.info.nodcount);
      NM = tmp_normals*mort_M; tmpmap = reshape([3*cont_face{2}.map.nod-2 3*cont_face{2}.map.nod-1 3*cont_face{2}.map.nod]',[],1); [i_,j_,v_] = find(NM); NM = sparse(i_,tmpmap(j_),v_,max(i_),mesh{2}.info.noddof*mesh{2}.info.nodcount);

      
      [new_activ_loc,new_inact_loc] = find_new_activ_inactiv(weight_gap, compl_param, z, normals_storage);
      
      sets_unchanged = isequal(new_activ_loc, activ_loc) && isequal(new_inact_loc, inact_loc);
      newton_continue = newton_continue || ~sets_unchanged;
      
      activ_loc = new_activ_loc; 
      inact_loc = new_inact_loc;
      inact_glob = cont_face{1}.map.nod(inact_loc);
      activ_glob = cont_face{1}.map.nod(activ_loc);
      
      k = k+1;
    end
  end
  export_ensight_result(     ensight_name, reshape(u{1}              ,[],1), struct(  'type','nod.vec'   ,'eletype',mesh{1}.ele.type,'name','Displacement'),sprintf('%03d',n)    );
  export_ensight_result_add( ensight_name, reshape(u{2}              ,[],1), struct(  'type','nod.vec'   ,'eletype',mesh{2}.ele.type,'name','Displacement'),sprintf('%03d',n),3,1);
  export_ensight_result(     ensight_name, reshape(bypr{1}.ele.sigma',[],1), struct(  'type','ele.tensym','eletype',mesh{1}.ele.type,'name','Stress'      ),sprintf('%03d',n)    );
  export_ensight_result_add( ensight_name, reshape(bypr{2}.ele.sigma',[],1), struct(  'type','ele.tensym','eletype',mesh{2}.ele.type,'name','Stress'      ),sprintf('%03d',n),3,1);


end