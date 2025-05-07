function plot_points_on_ref_element(element_type, points)
    figure;
    hold on;
    axis equal;
    grid on;
    xlabel('Xi'); ylabel('Eta');
    
    % Define reference elements
    switch element_type
        case 'tria3'
            ref_elem = [0, 0; 1, 0; 0, 1; 0, 0]; % Triangle + close the shape
        case 'quad4'
            ref_elem = [-1, -1; 1, -1; 1, 1; -1, 1; -1, -1]; % Quadrilateral + close
        otherwise
            error('Unsupported element type: %s', element_type);
    end
    
    % Plot the reference element
    plot(ref_elem(:,1), ref_elem(:,2), 'k-', 'LineWidth', 1.5);
    
    % Plot the points
    plot(points(:,1), points(:,2), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    
    title(sprintf('Reference Element: %s', element_type));
    hold off;
end



