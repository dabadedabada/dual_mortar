addpath('core_matfem');
addpath('core');
addpath('linearizations');

etype_s = 'hexa8';
etype_m = 'tetr4'; %'tetr4';
mesh_s = mesh_generator([6,6,1], [1,1,2],etype_s); 
mesh_m = mesh_generator([9,9,1], [3,3,2],etype_m);
mesh_s.nod.coo = mesh_s.nod.coo + [1.6,1.6,2];

cont_face_s = init_cont_face(mesh_s, 5);
cont_face_m = init_cont_face(mesh_m, 6);

FE = cell(2,1);
FE_s = fe_init(etype_s,2);
FE_m = fe_init(etype_m,2);

nn = cont_face_s.info.nN;
consts = cell(2,1);
consts_s = struct('Young', 1e3, 'Poisson', 0.3, 'density', 1e-1, 'lambda', [], 'mu', []);
consts_m = struct('Young', 2e3, 'Poisson', 0.4, 'density', 2e-1, 'lambda', [], 'mu', []);

consts_m.lambda = consts_m.Young*consts_m.Poisson/((1+consts_m.Poisson)*(1-2*consts_m.Poisson));   consts_m.mu     = consts_m.Young/(2*(1+consts_m.Poisson));
consts_s.lambda = consts_s.Young*consts_s.Poisson/((1+consts_s.Poisson)*(1-2*consts_s.Poisson));   consts_s.mu     = consts_s.Young/(2*(1+consts_s.Poisson));
compl_param = consts_s.Young;

normals_storage = get_normals_centroids(cont_face_s);
n0 = normals_storage.averaged_normals;
%n0       = normals_storage.ele_normals{1};                                 
dn0       = zeros([size(n0) size(n0,1) 3]);
n0_vla04 = vla04_getnormals(cont_face_s,zeros(3*size(cont_face_s.coo,1),1)); dn0_vla04 = dn0;
eeps = 1e-6;
for i = 1:size(n0,1)
  for j = 1:3
    tmp = cont_face_s;
    tmp.coo(i,j) = tmp.coo(i,j) + eeps;
    uij = zeros(3*nn,1); uij(3*(i-1)+j) = eeps;
    nij = get_normals_centroids(tmp);
    nij_vla04 = vla04_getnormals(cont_face_s,uij);
    dn0(      i,j,:,:) = (nij.averaged_normals-n0)/eeps;
    dn0_vla04(i,j,:,:) = (nij_vla04-n0_vla04)/eeps;
  end
end
Dn0       = zeros(nn,3,nn,3);
Dn0_vla04 = zeros(nn,3,nn,3);
for i = 1:nn
  tmp = linear_averaged_normals(cont_face_s, normals_storage, i);
  %   Dn0(:,:,1,1) = reshape(tmp(1,:),[nn, 3]);
  %   Dn0(:,:,1,2) = reshape(tmp(2,:),[nn, 3]);
  %   Dn0(:,:,1,3) = reshape(tmp(3,:),[nn, 3]);
  Dn0(:,:,i,1) = reshape(tmp(1,:),[3, nn])';
  Dn0(:,:,i,2) = reshape(tmp(2,:),[3, nn])';
  Dn0(:,:,i,3) = reshape(tmp(3,:),[3, nn])';
end
Dn0_vla04 = vla04_getDnormals(cont_face_s,zeros(3*size(cont_face_s.coo,1),1));


problem = cell(2,1);
problem{1} = struct('material_model', 'kirchhoff', 'volume_force', 0*[0,0,-10], 'variables', [], 'boundary', struct('dirichlet', []));
problem{1}.boundary.dirichlet = {...
  6,'0*[X,Y,Z]+[ 1.5,0,0]*max(T-0.2,0)+[0,0,-5.0]*min(T,0.2)','enforce'};
problem{2} = struct('material_model', 'kirchhoff', 'volume_force', 0*[0,0,-10], 'variables', [], 'boundary', struct('dirichlet', []));
problem{2}.boundary.dirichlet = {...
  1,'0*[X,Y,Z]+[ 0,0,0]','enforce'; ...
  2,'0*[X,Y,Z]+[ 0,0,0]','enforce'; ...
  3,'0*[X,Y,Z]+[ 0,0,0]','enforce'; ...
  4,'0*[X,Y,Z]+[ 0,0,0]','enforce'};
% problem{2}.boundary.dirichlet = {...
%   5,'0*[X,Y,Z]+[ 0,0,0]','enforce'};


figure(1);
hold on;
patch('Faces',[mesh_m.bou{1}.nod; mesh_m.bou{2}.nod; mesh_m.bou{3}.nod; ...
                mesh_m.bou{4}.nod; mesh_m.bou{5}.nod; mesh_m.bou{6}.nod],...
      'Vertices',mesh_m.nod.coo,'FaceColor','white');
patch('Faces',[mesh_s.bou{1}.nod; mesh_s.bou{2}.nod; mesh_s.bou{3}.nod; ...
                mesh_s.bou{4}.nod; mesh_s.bou{5}.nod; mesh_s.bou{6}.nod],...
      'Vertices',mesh_s.nod.coo,'FaceColor','white');


