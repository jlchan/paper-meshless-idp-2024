%paramters to circle

a = 3;
b = 3;
h = 0;
k = 0;
n_x = 25;
n_y = 25;
n_boundary = 75;
alpha = 2.5;
beta = 0.5;
scale_x = 1;
scale_y = 1;
epi = 1;

center_supp = zeros(1, 1);
n_boundary_supp = zeros(1, 1);

% %Uncomment below to use sharingan grid
% center_supp(1, 1) = 0;
% center_supp(1, 2) = 0;
% center_supp(1, 1) = 1.5*cos(2*pi/3);
% center_supp(1, 2) = 1.5*sin(2*pi/3);
% center_supp(2, 1) = 1.5*cos(4*pi/3); 
% center_supp(2, 2 )= 1.5*sin(4*pi/3);
% center_supp(3, 1) = 1.5*cos(6*pi/3);
% center_supp(3, 2) = 1.5*sin(6*pi/3);
% r_supp = zeros(1, 1);
% r_supp(1) = 2/3;
% r_supp(2) = 2/3;
% r_supp(3) = 2/3;

% n_boundary_supp = zeros(1, 1);
% n_boundary_supp(1) = 30;
% n_boundary_supp(2) = 30;
% n_boundary_supp(3) = 30;



%%
tic
options=optimoptions('quadprog','Algorithm','trust-region-reflective');
[points, H, H_min, Adj, L, E_x, E_y, S_x, S_y, Q_x, Q_y, phi_x, phi_y, normal_x, normal_y, elapsedTimeVec2] = matrix_generator(a, b, h, k, r_supp, center_supp, alpha, beta, n_boundary, n_boundary_supp, n_x, n_y, scale_x, scale_y);
elapsedTime = toc;

