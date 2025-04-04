function plot_points_3d(points)

% Plot the points in 3D
scatter3(points(:,1), points(:,2), points(:,3), 50, 'filled')

% Add labels and grid
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')
grid on
title('3D Scatter Plot')

end