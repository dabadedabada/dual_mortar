function plot_meshes(mesh)
% plots meshes created by mesh_generator

n = size(mesh, 1);
figure(1);
for i=1:n
  hold on;
  patch('Faces',[mesh{i}.bou{1}.nod; mesh{i}.bou{2}.nod; mesh{i}.bou{3}.nod; ...
                mesh{i}.bou{4}.nod; mesh{i}.bou{5}.nod; mesh{i}.bou{6}.nod],...
      'Vertices',mesh{i}.nod.coo,'FaceColor','white');
end
% axis properties for better visualization
axis equal;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);
hold off;