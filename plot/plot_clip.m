function plot_clip(x_s, y_s, x_m, y_m)
  % calculates the intersection between two polygons in 2D and plots it
  
  % Perform polygon intersection (METHOD = 1 for intersection)
  [X, Y] = polyclip(x_s, y_s, x_m, y_m, 1);
  
  % Plot the polygons and their intersection
  figure; hold on; axis equal;
  fill(x_s, y_s, 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 1.5); % Red polygon
  fill(x_m, y_m, 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 1.5); % Blue polygon
  % Plot intersection region (green)
  for i = 1:length(X)
      fill(X{i}, Y{i}, 'g', 'EdgeColor', 'k', 'LineWidth', 2);
  end
  
  legend('Slave', 'Master', 'Intersection');
  title('Polygon Clipping Example: Intersection');
  hold off;

end