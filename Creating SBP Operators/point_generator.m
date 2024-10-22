%Outputs points: A N x 3 matrix where N is the number of points in the
%pointi cloud. The first and second columns are the x and y coordinates
%respectively while the last column is binary which determines whetehr or
%not the point is on the boundary. The matrix is organized so that all the
%boundary points appear first before the interior points.

function points = point_generator(a, b, h, k, r_supp, center_supp, n_x, n_y, n_boundary, n_boundary_supp)
    left_end = h - a;
    lower_end = k - b;
    right_end = h + a;
    upper_end = k + b;
    N = n_x*n_y;
    
    
    points = zeros(N, 3);
    x_coords = linspace(left_end, right_end, n_x);
    y_coords = linspace(lower_end, upper_end, n_y);
    [x_coords, y_coords] = meshgrid(x_coords, y_coords);
    points(:, 1) = reshape(x_coords, [n_x*n_y, 1]);
    points(:, 2) = reshape(y_coords, [n_x*n_y, 1]);
    points(:, 3)  = in_ellipse(a, b, h, k, r_supp, center_supp, points);
    
    points = points(points(:, 3) ~= 0, :);
    
        
    
    if r_supp(1) > 0
        for i = length(r_supp):-1:1
            boundary_points_supp = zeros(n_boundary_supp(i)+1, 3);
            t_supp = linspace(0, 2*pi, n_boundary_supp(i)+1);
            x_coords_boundary_supp = r_supp(i)*cos(t_supp)+center_supp(i, 1);
            y_coords_boundary_supp = r_supp(i)*sin(t_supp)+center_supp(i, 2);
            boundary_points_supp(:, 1) = x_coords_boundary_supp;
            boundary_points_supp(:, 2) = y_coords_boundary_supp;
            boundary_points_supp(end, :) = [];
            points = [boundary_points_supp; points];

        end
    end
    
    boundary_points = zeros(n_boundary+1, 3);
    t = linspace(0, 2*pi, n_boundary+1);
    x_coords_boundary = a*cos(t)+h;
    y_coords_boundary = b*sin(t)+k;
    boundary_points(:, 1) = x_coords_boundary;
    boundary_points(:, 2) = y_coords_boundary;
    boundary_points(end, :) = [];
    points = [boundary_points ; points];


    
end

function [val] = in_ellipse(a, b, h, k, r_supp, center_supp, coord)
   if r_supp(1) > 0
       val_all = zeros(size(coord, 1), length(r_supp) + 1);
       val_all(:, 1) = (coord(:, 1) - h).^2 / a^2 + (coord(:, 2) - k).^2 / b^2 < 1;
       for i = 1:length(r_supp)
           val_all(:, i+1) = (coord(:, 1) - center_supp(i, 1)).^2  + (coord(:, 2) - center_supp(i, 2)).^2  > r_supp(i)^2;
       end

       val = all([val_all], 2);
   else
       val = (coord(:, 1) - h).^2 / a^2 + (coord(:, 2) - k).^2 / b^2 < 1;
   end

end