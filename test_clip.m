clear
addpath(genpath('..'));

% Define first polygon (Quadrilateral)
x1 = [0  2  2.5  0.5  0];  % X-coordinates
y1 = [0  0  2  1.5  0];  % Y-coordinates

% Define second polygon (Quadrilateral, overlapping)
x2 = [1  3  exp(1)  1  1];  % X-coordinates
y2 = [1  1  exp(1)  3  1];  % Y-coordinates

% Perform polygon intersection (METHOD = 1 for intersection)
[X, Y] = polyclip(x1, y1, x2, y2, 1);

% Plot the polygons and their intersection
figure; hold on; axis equal;
fill(x1, y1, 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 1.5); % Red polygon
fill(x2, y2, 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 1.5); % Blue polygon

% Plot intersection region (green)
for i = 1:length(X)
    fill(X{i}, Y{i}, 'g', 'EdgeColor', 'k', 'LineWidth', 2);
end

legend('Polygon 1', 'Polygon 2', 'Intersection');
title('Polygon Clipping Example: Intersection');
hold off;
