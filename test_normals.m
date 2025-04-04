clear
addpath(genpath('..'));

% Create reference configurations X^(1,2), t=0
mesh_s = mesh_generator([6,6,1], [3,3,2],'hexa8'); 
mesh_m = mesh_generator([9,9,1], [3,3,2],'tetr4');
mesh_s.nod.coo=mesh_s.nod.coo + [0,0,5];
mesh_s.nod.coo=mesh_s.nod.coo + [0,2,0];
mesh_s.nod.coo=mesh_s.nod.coo + [2,0,0];

% Changes to the geometry
mesh_s.nod.coo(:, 3) = mesh_s.nod.coo(:, 3) - mesh_s.nod.coo(:, 1) / 5;
mesh_s.nod.coo(mesh_s.bou{5}.nod(1,1), :) = mesh_s.nod.coo(mesh_s.bou{5}.nod(1,1), :) - 0.5;
mesh_s.nod.coo(mesh_s.bou{5}.nod(1,2), :) = mesh_s.nod.coo(mesh_s.bou{5}.nod(1,2), :) - 0.5;
%mesh_m.nod.coo(mesh_m.bou{6}.nod(1,2), :) = mesh_m.nod.coo(mesh_m.bou{6}.nod(1,2), :) - 0.5;


figure(1);
hold on;
patch('Faces',[mesh_m.bou{1}.nod; mesh_m.bou{2}.nod; mesh_m.bou{3}.nod; ...
                mesh_m.bou{4}.nod; mesh_m.bou{5}.nod; mesh_m.bou{6}.nod],...
      'Vertices',mesh_m.nod.coo,'FaceColor','white');
patch('Faces',[mesh_s.bou{1}.nod; mesh_s.bou{2}.nod; mesh_s.bou{3}.nod; ...
                mesh_s.bou{4}.nod; mesh_s.bou{5}.nod; mesh_s.bou{6}.nod],...
      'Vertices',mesh_s.nod.coo,'FaceColor','white');

% axis properties for better visualization
axis equal;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);
hold off;


% init contact faces
cont_face_s = init_cont_face(mesh_s, 5);
cont_face_m = init_cont_face(mesh_m, 6);

% Sizes
nele_s = cont_face_s.info.nele;
nele_m = cont_face_m.info.nele;
nN_s_ele = cont_face_s.info.nN_ele;
nN_m_ele = cont_face_m.info.nN_ele;
nN_s = cont_face_s.info.nN;
nN_m = cont_face_m.info.nN;

% Get the normals by sum of nodal normals and then mapping at the center
% not really averaged, just normalized sum

% Normals
normals_storage = get_normals_centroids(cont_face_s);
averaged_normals = zeros(size(normals_storage.sum_normals));
for i=1:nN_s
  averaged_normals(i,:) = normals_storage.sum_normals(i,:)/norm(normals_storage.sum_normals(i,:));
end

% Visualize nodal averaged normals and tangents
figure(1);
hold on;
quiver3(cont_face_s.coo(:,1)', cont_face_s.coo(:,2)', cont_face_s.coo(:,3)', averaged_normals(:,1)', averaged_normals(:,2)', averaged_normals(:,3)', 1, 'r', 'LineWidth', 1);
quiver3(cont_face_s.coo(:,1)', cont_face_s.coo(:,2)', cont_face_s.coo(:,3)', normals_storage.t1(:,1)', normals_storage.t1(:,2)', normals_storage.t1(:,3)', 1, 'r', 'LineWidth', 1);
quiver3(cont_face_s.coo(:,1)', cont_face_s.coo(:,2)', cont_face_s.coo(:,3)', normals_storage.t2(:,1)', normals_storage.t2(:,2)', normals_storage.t2(:,3)', 1, 'r', 'LineWidth', 1);
legend({'Slave Mesh', 'Slave Normals'});
hold off;







