%Function point_generator creates a N x 3 matrix where N is the total number
%of points. The first column is the x coordintes, the 2nd column is the y 
%coordinates and the last column is binary that determines whether the
%point is on the boudnary.

function points = point_generator(a, b, h, k, n_x, n_y, n_boundary)
    left_end = h - a;
    lower_end = k - b;
    right_end = h + a;
    upper_end = k + b;
    N = n_x*n_y;

    
    
    points = zeros(N, 4);
    x_coords = linspace(left_end, right_end, n_x);
    y_coords = linspace(lower_end, upper_end, n_y);
    [x_coords, y_coords] = meshgrid(x_coords, y_coords);
    points(:, 1) = reshape(x_coords, [n_x*n_y, 1]);
    points(:, 2) = reshape(y_coords, [n_x*n_y, 1]);
    points(:, 3)  = in_ellipse(a, b, h, k, points);
    
    points = points(points(:, 3) ~= 0, :);
    
    
    
    points(:, 4) = [n_boundary + 1:n_boundary + size(points, 1)];
    
    boundary_points = zeros(n_boundary+1, 4);
    t = linspace(0, 2*pi, n_boundary+1);
    x_coords_boundary = a*cos(t)+h;
    y_coords_boundary = b*sin(t)+k;
    boundary_points(:, 1) = x_coords_boundary;
    boundary_points(:, 2) = y_coords_boundary;
    boundary_points(:, 4) = [1:n_boundary+1];
    boundary_points(end, :) = [];
    points = [boundary_points ; points];
    
    
    
end

function [val] = in_ellipse(a, b, h, k, coord)
   val = (coord(:, 1) - h).^2 / a^2 + (coord(:, 2) - k).^2 / b^2 < 1;
   
end
