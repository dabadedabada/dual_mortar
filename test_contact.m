clear
%addpath(genpath('..'));
addpath("plot");
addpath('core_matfem');
addpath('core');
addpath('clipper2');
%addpath('../../../../programming/matlab/matfem_marta/core');
addpath("linearizations");

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
problem{1} = struct('material_model', 'kirchhoff', 'volume_force', 0*[0,0,-10], 'variables', [], 'boundary', struct('dirichlet', []));
problem{1}.boundary.dirichlet = {...
  6,'0*[X,Y,Z]+[ 1.5,0,0]*max(T-0.2,0)+[0,0,-3.0]*min(T,0.2)','enforce'};
problem{2} = struct('material_model', 'kirchhoff', 'volume_force', 0*[0,0,-10], 'variables', [], 'boundary', struct('dirichlet', []));
%problem{2}.boundary.dirichlet = {...
%  1,'0*[X,Y,Z]+[ 0,0,0]','enforce'; ...
%  2,'0*[X,Y,Z]+[ 0,0,0]','enforce'; ...
%  3,'0*[X,Y,Z]+[ 0,0,0]','enforce'; ...
%  4,'0*[X,Y,Z]+[ 0,0,0]','enforce'};
problem{2}.boundary.dirichlet = {...
  5,'0*[X,Y,Z]+[ 0,0,0]','enforce'};

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


inact_loc = 1:size(cont_face{1}.map.nod, 1);       % init inactive 0 2D indices
activ_loc = [];                                    % init active 0   2D indices

inact_glob = cont_face{1}.map.nod(inact_loc)';      % init inactive 0 3D indices
activ_glob = cont_face{1}.map.nod(activ_loc)';      % init active 0   3D indices

m_id_loc = 1:size(cont_face{2}.map.nod, 1);        % init set{M} 2D indices
idx_m_loc = reshape([3*m_id_loc-2; 3*m_id_loc-1; 3*m_id_loc], 1, []);
m_id_glob = cont_face{2}.map.nod';                  % init set{M} 3D indices
idx_m_glob = reshape([3*m_id_glob-2; 3*m_id_glob-1; 3*m_id_glob], 1, []) + 3*mesh{1}.info.nodcount;

s_id_glob = cont_face{1}.map.nod';
idx_s_glob = reshape([3*s_id_glob-2; 3*s_id_glob-1; 3*s_id_glob], 1, []);

z = zeros(3*cont_face{1}.info.nN,1);       % init z0

% ---------- Reference configuration -------------
% Normals, centroids and tangents
normals_storage = get_normals_centroids(cont_face{1});

% Algorithm 1 - for reference configuration X
[mort_D, mort_M, weight_gap, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cont_face{1}, cont_face{2}, normals_storage);
[mort_Dalg, mort_Malg, weight_gap_alg] = linear_dual_mortar_fem(cont_face{1}, cont_face{2}, clips_storage, normals_storage, slave_storage, master_storage);
%weight_gap_dif = weight_gap_alg;
tmp_normals = mat2cell(normals_storage.averaged_normals,ones(size(normals_storage.averaged_normals,1),1),3);
tmp_normals = blkdiag(tmp_normals{:});

% ---------- time scheme -----------------------

for n = 1:timestepping.steps
  t  =   timestepping.times(n+1);
  dt = t-timestepping.times(n);

  fprintf('Time (%3i) %f \n',n,t);
  BC{1} = set_boundary_conditions(mesh{1}, problem{1}, t, n, u{1} );
  BC{2} = set_boundary_conditions(mesh{2}, problem{2}, t, n, u{2} );

  k = 1;
  newton_continue = true;
  while newton_continue
    [disc{1}, bypr{1}] = assemble_hyperelasticity(FE{1}, mesh{1}, consts{1}, problem{1}, u{1});
    [disc{2}, bypr{2}] = assemble_hyperelasticity(FE{2}, mesh{2}, consts{2}, problem{2}, u{2});

    %if n==1 && k==1
    %  K = blkdiag(disc{1}.Kc, disc{2}.Kc);
    %end
    f_c = zeros(3*(mesh{1}.info.nodcount + mesh{2}.info.nodcount),1);
    f_c(idx_s_glob) = mort_D*z; f_c(idx_m_glob) = -mort_M'*z;

    DK = blkdiag( disc{1}.Kc+disc{1}.Ks+(0*dt^-2)*disc{1}.M,                disc{2}.Kc+disc{2}.Ks    +(0*dt^-2)*disc{2}.M);

    r  =       [(-disc{1}.fext+disc{1}.fint   +(0*dt^-2)*disc{1}.M*u{1}); (-disc{2}.fext+disc{2}.fint+(0*dt^-2)*disc{2}.M*u{2})] + f_c;
    %DK = K;
    %r  =       [-disc{1}.fext; -disc{2}.fext] + f_c + K*[u{1};u{2}];
    idx_inact_loc = reshape([3*inact_loc-2; 3*inact_loc-1; 3*inact_loc], 1, []);
    idx_activ_loc = reshape([3*activ_loc-2; 3*activ_loc-1; 3*activ_loc], 1, []);

    idx_inact_glob = reshape([3*inact_glob-2; 3*inact_glob-1; 3*inact_glob], 1, []);
    idx_activ_glob = reshape([3*activ_glob-2; 3*activ_glob-1; 3*activ_glob], 1, []);

    Ctilde = [squeeze(pagemtimes(permute(-mort_Malg, [2, 1, 3]),z));squeeze(pagemtimes(mort_Dalg,z))];
    Ctilde_MM = Ctilde(idx_m_loc,idx_m_loc);
    Ctilde_II = Ctilde(idx_inact_loc+ 3*cont_face{2}.info.nN,idx_inact_loc+ 3*cont_face{2}.info.nN);
    Ctilde_AA = Ctilde(idx_activ_loc+ 3*cont_face{2}.info.nN,idx_activ_loc+ 3*cont_face{2}.info.nN);

    Ctilde_AM = Ctilde(idx_activ_loc+ 3*cont_face{2}.info.nN,idx_m_loc);
    Ctilde_AI = Ctilde(idx_activ_loc+ 3*cont_face{2}.info.nN,idx_inact_loc+ 3*cont_face{2}.info.nN);
    Ctilde_IM = Ctilde(idx_inact_loc+ 3*cont_face{2}.info.nN,idx_m_loc);

    DK(idx_m_glob ,idx_m_glob ) = DK(idx_m_glob ,idx_m_glob ) + Ctilde_MM;
    DK(idx_activ_glob,idx_activ_glob) = DK(idx_activ_glob,idx_activ_glob) + Ctilde_AA;
    DK(idx_inact_glob,idx_inact_glob) = DK(idx_inact_glob,idx_inact_glob) + Ctilde_II;
    DK(idx_activ_glob,idx_m_glob ) = DK(idx_activ_glob,idx_m_glob ) + Ctilde_AM;   DK(idx_m_glob ,idx_activ_glob) = DK(idx_m_glob ,idx_activ_glob) + Ctilde_AM';
    DK(idx_inact_glob,idx_m_glob ) = DK(idx_inact_glob,idx_m_glob ) + Ctilde_IM;   DK(idx_m_glob ,idx_inact_glob) = DK(idx_m_glob ,idx_inact_glob) + Ctilde_IM';
    DK(idx_activ_glob,idx_inact_glob) = DK(idx_activ_glob,idx_inact_glob) + Ctilde_AI;   DK(idx_inact_glob,idx_activ_glob) = DK(idx_inact_glob,idx_activ_glob) + Ctilde_AI';

    Z1 = zeros(3*(mesh{1}.info.nodcount + mesh{2}.info.nodcount), 3*cont_face{1}.info.nN);
    Z1(idx_s_glob,:) = mort_D; Z1(idx_m_glob,:) = -mort_M';

    row5 = zeros(length(idx_inact_loc), 3*(mesh{1}.info.nodcount + mesh{2}.info.nodcount));
    Z2 = zeros(length(idx_inact_loc), 3*cont_face{1}.info.nN);
    Z2(:,idx_inact_loc) = eye(length(idx_inact_loc));

    row6 = zeros(length(activ_loc), 3*(mesh{1}.info.nodcount + mesh{2}.info.nodcount));
    Stilde_A = weight_gap_alg(activ_loc, 3*cont_face{2}.info.nN+1:end);
    Mtilde_A =  weight_gap_alg(activ_loc, 1:3*cont_face{2}.info.nN);

    %Stilde_A = weight_gap_dif(activ_loc, 3*cont_face{2}.info.nN+1:end);
    %Mtilde_A = weight_gap_dif(activ_loc, 1:3*cont_face{2}.info.nN);
    row6(:, idx_m_glob) = Mtilde_A; row6(:, idx_s_glob) = Stilde_A;


    row7 = zeros(2*length(activ_loc),3*(mesh{1}.info.nodcount + mesh{2}.info.nodcount));
    Z3 = zeros(2*length(activ_loc), 3*cont_face{1}.info.nN);
    if ~isempty(activ_loc)
      zi = reshape(z, 3, [])';
      for i=1:length(activ_loc)
        [T1, T2] = linear_get_tangents(cont_face{1}, normals_storage, activ_loc(i));
        row7(i,idx_s_glob) = zi(activ_loc(i),:)*T1;
        row7(length(activ_loc)+i,idx_s_glob) = zi(activ_loc(i),:)*T2;
        Z3(i,3*activ_loc(i)-2:3*activ_loc(i)) = normals_storage.t1(activ_loc(i),:);
        Z3(length(activ_loc)+i,3*activ_loc(i)-2:3*activ_loc(i)) = normals_storage.t2(activ_loc(i),:);
      end
    end
    % solve the system
    B_ = [blkdiag(BC{1}.Bd{1},BC{2}.Bd{1})];
    c_ = [BC{1}.Bd{1}*(BC{1}.U0-u{1}); BC{2}.Bd{1}*(BC{2}.U0-u{2})];

    % condense dirichlet to DK and r
    nB = size(B_,1);
    [~,indEQU] = find(B_); logindEQU = false(nB,1); logindEQU(indEQU) = true;
    sc_  = 1;
    DK_  = DK/sc_;
    r_   = r/sc_;
    DK_d = DK_(~logindEQU,logindEQU)*c_;
    DK_(:,logindEQU)=0; DK_(logindEQU,:)=0; DK_(logindEQU,logindEQU)=speye(nB);
    r_( logindEQU,:)=-c_; r_(~logindEQU,:)=r_(~logindEQU,:)+DK_d;

    AA = [DK_ Z1; ...
      row5 Z2; ...
      row6 zeros(size(row6,1),size(Z2,2)); ...
      row7 Z3
      ];
    bb = [-r_; ...
      zeros(size(row5,1),1); ...
      -weight_gap(activ_loc); ...
      zeros(size(row7,1),1)...
      ];
    result = AA\bb;

    du{1} = result(             1:size(u{1},1));
    du{2} = result(size(u{1},1)+1:size(u{1},1)+size(u{2},1));
    z = result(end-3*cont_face{1}.info.nN+1:end);
    if k == 1,      newton_err = norm([du{1}; du{2}]);
    else,           newton_err = norm([du{1}; du{2}])/norm([u{1}; u{2}]);
    end
    fprintf('  newton %3d: %1.8e  [A/I]=[%d/%d] \n',k,newton_err, length(activ_glob), length(inact_glob) );
    newton_continue = (k <= 1) || (newton_err > 1e-6);

    % test KKT
    %z_inact_should_be_zero = [row5 Z2]*result;                             % should be zero
    %wgap_equals_its_linear = [row6 zeros(size(row6,1),size(Z2,2))]*result; % = -weight_gap(activ_loc)
    %t1t2_should_cancelout = [row7 Z3]*result;                              % should be almost zero

    % test DK
    %testfc = [DK_ Z1]*result; % -r_(idx_m_glob) == fc(idx_m_glob)

    % update displacement
    u{1} = u{1} + du{1};
    u{2} = u{2} + du{2};
    
    % add d^k+1 to cont_face coordinates
    cont_face{1}.coo = mesh{1}.nod.coo(cont_face{1}.map.nod,:) + reshape(u{1}(idx_s_glob),3,[])';
    cont_face{2}.coo = mesh{2}.nod.coo(cont_face{2}.map.nod,:) + reshape(u{2}(idx_m_glob-3*mesh{1}.info.nodcount),3,[])';

    % Normals, centroids and tangents for d^k+1
    normals_storage = get_normals_centroids(cont_face{1});

    % Algorithm 1 - for current configuration with d^k+1
    [mort_D, mort_M, weight_gap, clips_storage, slave_storage, master_storage] = get_contact_dual_mortar(cont_face{1}, cont_face{2}, normals_storage);
    [mort_Dalg, mort_Malg, weight_gap_alg] = linear_dual_mortar_fem(cont_face{1}, cont_face{2}, clips_storage, normals_storage, slave_storage, master_storage);

    [new_activ_loc,new_inact_loc] = find_new_activ_inactiv(weight_gap, compl_param, z, normals_storage);

    sets_unchanged = isequal(new_activ_loc, activ_loc) && isequal(new_inact_loc, inact_loc);
    %newton_continue =  ~(~newton_continue && ~sets_unchanged);

    activ_loc = new_activ_loc';
    inact_loc = new_inact_loc';
    inact_glob = cont_face{1}.map.nod(inact_loc)';
    activ_glob = cont_face{1}.map.nod(activ_loc)';

    k = k+1;
  end
  export_ensight_result(     ensight_name, reshape(u{1}              ,[],1), struct(  'type','nod.vec'   ,'eletype',mesh{1}.ele.type,'name','Displacement'),sprintf('%03d',n)    );
  export_ensight_result_add( ensight_name, reshape(u{2}              ,[],1), struct(  'type','nod.vec'   ,'eletype',mesh{2}.ele.type,'name','Displacement'),sprintf('%03d',n),3,1);
  export_ensight_result(     ensight_name, reshape(bypr{1}.ele.sigma',[],1), struct(  'type','ele.tensym','eletype',mesh{1}.ele.type,'name','Stress'      ),sprintf('%03d',n)    );
  export_ensight_result_add( ensight_name, reshape(bypr{2}.ele.sigma',[],1), struct(  'type','ele.tensym','eletype',mesh{2}.ele.type,'name','Stress'      ),sprintf('%03d',n),3,1);
end

