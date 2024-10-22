%Function outpus the various matrices needed to create the SBP operators
%H: the uniform norm matrix
%H_min: the optimized H matrix
%A, L: adjacency and graph laplacian matrix
%E_x, E_y, S_x, S_y, Q_x, Q_y: the SBP matrices
%phi_x, phi_y: vectors used to create S_x and S_y
%normal_x, normal_y: normal vectors at the boundary



function [points, H, H_min, A, L, E_x, E_y, S_x, S_y, Q_x, Q_y, phi_x, phi_y, normal_x, normal_y, elapsedTimeVec] = matrix_generator_optimized(a, b, h, k, r_supp, center_supp, alpha, beta, n_boundary,n_boundary_supp, n_x, n_y, scale_x, scale_y) 
    elapsedTimeVec = [];
    tic;
    points = point_generator(a, b, h, k, r_supp, center_supp, n_x, n_y, n_boundary, n_boundary_supp);
    elapsedTimeVec(1) = toc;

    tic;
    H = H_matrix(points, a, b);
    elapsedTimeVec(2) = toc;
    
    tic;
    %minspan adjacency
    %A = adjacency_matrix_min_span(points, D, a, b, h, k, alpha, beta, n_x, n_y);
    
    %degree 1/degree 2 delauntay triangulation adjacency
    %A = adjacency_matrix_triangulation(points, a, b, h, k, alpha, beta, n_x, n_y, n_boundary, n_boundary_supp);
    %elapsedTimeVec(4) = toc;
    
    %Euclidean radius adjacency
    A = adjacency_matrix_KT(points, a, b, h, k, alpha, beta, n_x, n_y, n_boundary, n_boundary_supp, scale_x, scale_y);
    elapsedTimeVec(4) = toc;

    tic;
    L = laplacian_matrix(points, A);
    elapsedTimeVec(5) = toc;
    
    tic;
    [E_x, E_y, normal_x, normal_y] = E_matrix(points, a, b, h, k, n_boundary, r_supp, center_supp, n_boundary_supp);
    elapsedTimeVec(6) = toc;

    tic;
    [phi_x, phi_y] = potential_vector(points, L, E_x, E_y);
    elapsedTimeVec(7) = toc;

    tic;
    [S_x, S_y] = S_matrix(points, phi_x, phi_y, L);
    elapsedTimeVec(8) = toc;

    tic;
    Q_x = S_x + 0.5*E_x;
    Q_y = S_y + 0.5*E_y;
    elapsedTimeVec(9) = toc;
    
    tic;
    H_min = H_minimizer(points, Q_x, Q_y, H);
    elapsedTimeVec(10) = toc;
end

function H = H_matrix(points, a, b)
    diagonal_elements = zeros(size(points, 1), 1) + (pi*a*b / size(points, 1));
    H = spdiags(diagonal_elements, 0, length(diagonal_elements), length(diagonal_elements));



end

function H_min = H_minimizer(points, Q_x, Q_y, H)
   addpath("/Users/SamuelKwan/Documents/Rice/Research CAAM/Meshless/scs-matlab")
   N = size(points, 1);
   x_coords = points(:, 1);
   y_coords = points(:, 2);
   D = sparse(2*speye(N));
   f = full(-(Q_x*x_coords + Q_y*y_coords));
   epsilon = 1/N^2;
   lb = zeros(N, 1) + epsilon;
   
   data.P = D;
   data.c = f;
   data.A = -speye(N);
   data.b = -lb;
   data.x = full(diag(H));
   cone.l = N;
   
   % % Optional solver settings
   settings = struct('eps_abs', 1e-13, 'eps_rel', 1e-13, 'eps_infeas', 1e-13);
   
   
   [H_min, ~, ~, ~] = scs(data, cone, settings);

   
end 


function A = adjacency_matrix_min_span(points, dist_matrix, a, b, h, k, alpha, beta, n_x, n_y)
    N = size(points, 1);
    A = zeros(N);
    G = graph(dist_matrix);
    min_span_tree = minspantree(G);
    min_span_tree_indexed = table2array(min_span_tree.Edges);
    
    for i = 1:size(min_span_tree_indexed, 1)
        A(min_span_tree_indexed(i, 1), min_span_tree_indexed(i, 2)) = 1;
    end
    A = A.' + A;

    
end

function A = adjacency_matrix_triangulation(points, a, b, h, k, alpha, beta, n_x, n_y, n_boundary, n_boundary_supp)
    N = size(points, 1);
    DT = delaunayTriangulation(points(:, 1:2));
    Edges = DT.edges;
    A = sparse(Edges(:,1), Edges(:,2), 1, N, N);
    A = A.' + A;
    
   %%Uncomment below for degree 2
   %A(1:n_boundary + sum(n_boundary_supp), 1:n_boundary + sum(n_boundary_supp)) = 0;
   %A = A + A^2;
    A = A - diag(diag(A));
    A = A > 0;
    A(1:n_boundary + sum(n_boundary_supp), 1:n_boundary + sum(n_boundary_supp)) = 0;
    
 
end

function A = adjacency_matrix_KT(points, a, b, h, k, alpha, beta, n_x, n_y, n_boundary, n_boundary_supp, scale_x, scale_y)
    N = size(points, 1);
    Mdl = KDTreeSearcher(points(:, 1:2), 'Distance','Euclidean');
    dx = 2*a / (n_x - 1);
    dy = 2*b / (n_y - 1);
    Idx = rangesearch(Mdl, points(:, 1:2), alpha*max(scale_x*dx, scale_y*dy), 'Distance','Euclidean');
    vectorLengths = cellfun(@numel, Idx);
    A = sparse(repelem(1:N, vectorLengths).' , cell2mat(Idx.'), 1, N, N);
    A = A - diag(diag(A));

end

function L = laplacian_matrix(points, A)
    N = size(points, 1);
    G = graph(A);
    L = sparse(laplacian(G));
    
end

function [phi_x, phi_y] = potential_vector(points, L, E_x, E_y) 
    N = size(points, 1);
    A_mat_1 = sparse(horzcat(vertcat(L, ones(1, N)), [ones(N, 1); 0]));
    A_mat_2 = sparse(horzcat(vertcat(L, ones(1, N)), [ones(N, 1); 0]));
    b_vec_1 = sparse(0.5*[E_x*ones(N, 1); 0]);
    b_vec_2 = sparse(0.5*[E_y*ones(N, 1); 0]); 
    
    phi_x_tmp = A_mat_1\b_vec_1;
    phi_y_tmp = A_mat_2\b_vec_2;
    phi_x = phi_x_tmp(1:end-1);
    phi_y = phi_y_tmp(1:end-1);

        
    
end


function [E_x, E_y, normal_x, normal_y] = E_matrix(points, a, b, h, k, n_boundary, r_supp, center_supp, n_boundary_supp)
    t = linspace(0, 2*pi, n_boundary+1);
    delta_t = t(2) - t(1);
    w = r_prime_norm(a, b, t)*delta_t;
    w(:, end) = [];
    [normal_x, normal_y] = normal_calc(points, a, b, h, k, 1, n_boundary);
    
    if r_supp(1) > 0
        start = n_boundary + 1;
        for i = 1:length(r_supp)
            if i > 1
                start = start + n_boundary_supp(i-1);
            end
            t = linspace(0, 2*pi, n_boundary_supp(i)+1);
            delta_t = t(2) - t(1);
            w_supp = r_prime_norm(r_supp(i), r_supp(i), t)*delta_t;
            w_supp(:, end) = [];
            [normal_x_supp, normal_y_supp] = normal_calc(points, r_supp(i), r_supp(i), center_supp(i, 1), center_supp(i, 2), start, n_boundary_supp(i));
            normal_x = [normal_x; -normal_x_supp];
            normal_y = [normal_y; -normal_y_supp];
            w = [w, w_supp];

        end
    end
     N = size(points, 1);



    idx = 1:n_boundary + sum(n_boundary_supp);
    
    
    E_x = sparse(idx, idx, normal_x .* w.', N, N);
    E_y = sparse(idx, idx, normal_y .* w.', N, N);
    
    
    
end

function val = r_prime_norm(a, b, t)
    val = sqrt(a^2*sin(t).^2+b^2*cos(t).^2);
end

function [normal_x, normal_y] = normal_calc(points, a, b, h, k, start, n_boundary)
     normal_y = a^2*(points([start:start + n_boundary - 1], 2) - k);
     normal_x = b^2*(points([start:start + n_boundary - 1], 1) - h);
     norm = sqrt(normal_y.^2 + normal_x.^2);
     normal_y = normal_y ./ norm;
     normal_x = normal_x ./ norm;
    
end

function [S_x, S_y] = S_matrix(points, phi_x, phi_y, L) 
    N = size(points, 1);
   [row, col] = find(abs(L) > 1e-12);
    


    S_x = sparse(row, col, phi_x(col) - phi_x(row));
    S_y = sparse(row, col, phi_y(col) - phi_y(row));


end  