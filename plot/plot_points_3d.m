function plot_points_3d(points)
% Plot the points in 3D with indices

% Scatter plot
scatter3(points(:,1), points(:,2), points(:,3), 50, 'filled')
hold on

% Label each point with its index
for i = 1:size(points,1)
    text(points(i,1), points(i,2), points(i,3), sprintf('  %d', i), ...
        'FontSize', 10, 'Color', 'white')
end

% Add labels and grid
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')
grid on
title('3D Scatter Plot')

hold off
end