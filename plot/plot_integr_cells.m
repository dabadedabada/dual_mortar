function plot_integr_cells(integr_cells)
  ncells = length(integr_cells); 

  figure;
  hold on;
  axis equal;
  xlabel('X'); ylabel('Y'); zlabel('Z');
  title('Triangular Integration Cells');
  
  for i = 1:ncells
      % Extract the 3x3 matrix for the i-th cell
      tri = integr_cells{i}; 
      
      % Plot the triangle using patch
      patch('Vertices', tri, 'Faces', [1 2 3], 'FaceColor', 'cyan', 'EdgeColor', 'black', 'FaceAlpha', 0.5);
  end
  
  hold off;
  view(3); % Set 3D view
end