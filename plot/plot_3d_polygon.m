function plot_3d_polygon(points1, points2, points3)
    % PLOT_3D_POLYGON Visualizes a 3D polygon using patch()
    %
    % INPUTS:
    %   points - Nx3 matrix of 3D coordinates [x, y, z]
    %   faces  - MxK matrix, each row contains indices of 'points' defining a face
    %
    % EXAMPLE USAGE:
    %   points = [0 0 0; 1 0 0; 1 1 1; 0 1 1];
    %   faces = [1 2 3 4];
    %   plot_3d_polygon(points, faces);

    faces1 = 1:size(points1, 1);
    faces2 = 1:size(points2, 1);
    faces3 = 1:size(points3, 1);

    figure; hold on; grid on; axis equal;
    
    % Create the polygon using patch
    patch('Vertices', points1, 'Faces', faces1, 'FaceColor', 'cyan', 'FaceAlpha', 0.5);
    patch('Vertices', points2, 'Faces', faces2, 'FaceColor', 'red', 'FaceAlpha', 0.5);
    patch('Vertices', points3, 'Faces', faces3, 'FaceColor', 'yellow', 'FaceAlpha', 0.5);
    
    % Mark the vertices with red circles
    scatter3(points1(:,1), points1(:,2), points1(:,3), 'ro', 'filled');
    scatter3(points2(:,1), points2(:,2), points2(:,3), 'ro', 'filled');
    scatter3(points3(:,1), points3(:,2), points3(:,3), 'ro', 'filled');
    
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('3D Polygon Visualization');
    view(3);
    
    hold off;
end