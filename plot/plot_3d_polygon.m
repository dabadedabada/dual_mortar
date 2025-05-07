function plot_3d_polygon(points1, points2)


    faces1 = 1:size(points1, 1);
    faces2 = 1:size(points2, 1);

    figure; hold on; grid on; axis equal;
    
    % Create the polygon using patch
    patch('Vertices', points1, 'Faces', faces1, 'FaceColor', 'cyan', 'FaceAlpha', 0.5);
    patch('Vertices', points2, 'Faces', faces2, 'FaceColor', 'red', 'FaceAlpha', 0.5);
    
    % Mark the vertices with red circles
    scatter3(points1(:,1), points1(:,2), points1(:,3), 'ro', 'filled');
    scatter3(points2(:,1), points2(:,2), points2(:,3), 'ro', 'filled');
    
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('3D Polygon Visualization');
    view(3);
    
    hold off;
end