%%
%paramters to circle

a = 3;
b = 3;
h = 0;
k = 0;
n_x = 25;
n_y = 25;
n_boundary = 75;
alpha = 2.5;
scale_x = 1;
scale_y = 1;
epi = 1;

center_supp = zeros(1, 1);
n_boundary_supp = zeros(1, 1);

% %Uncomment below to use sharingan point cloud
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
%creating the operators in case it isnt done so yet
options=optimoptions('quadprog','Algorithm','trust-region-reflective');
[points, H, H_min, Adj, L, E_x, E_y, S_x, S_y, Q_x, Q_y, phi_x, phi_y, normal_x, normal_y, elapsedTimeVec2] = matrix_generator(a, b, h, k, r_supp, center_supp, alpha, n_boundary, n_boundary_supp, n_x, n_y, scale_x, scale_y);
%%
%functions used are f(x, y) = 1 and f(x, y) = 4*sin(x)*sin(y)*e^(-(x^2+y^2))
N = size(points, 1);


%testing f(x, y) = 1
x_diff_1 = ((1./H_min).*(Q_x*points(:, 1)))-1;
y_diff_1 = ((1./H_min).*(Q_y*points(:, 2))) - 1;


%testing f(x, y) = 4*sin(x)*sin(y)*e^(-(x^2+y^2))
x_diff_func = (1./H_min).*(Q_x*func(points(:, 1), points(:, 2))) - dx_func(points(:, 1), points(:, 2));
y_diff_func = (1./H_min).*(Q_y*func(points(:, 1), points(:, 2))) - dy_func(points(:, 1), points(:, 2));


%%
%Displays the L2 error

sqrt(sum(H_min.*(x_diff_1.^2)))
sqrt(sum(H_min.*(y_diff_1.^2)))
sqrt(sum(H_min.*(x_diff_func.^2)))
sqrt(sum(H_min.*(y_diff_func.^2)))

%%
%Interior point vs exterior point experiment.
n_x = 150;
n_b = 500;
n_i = 0;

ext = find(points(:, 1).^2 + points(:, 2).^2 >= 4);
int = find(points(:, 1).^2 + points(:, 2).^2 < 4);


sqrt(max(abs(x_diff_1(ext).^2)))
sqrt(max(abs(y_diff_1(ext).^2)))

sqrt(max(abs(x_diff_1(int).^2)))
sqrt(max(abs(y_diff_1(int).^2)))


%%
%plotting the error heat maps


figure()
scatter(points(:, 1), points(:, 2), [], abs(x_diff_func), 'filled');
colormap(gca, 'parula');
cb = colorbar;
set(cb, 'FontSize', 35);
cb.Position = [0.760, 0.1, 0.032 0.8]; % [left, bottom, width, height]
ax = gca; % Get current axes
ax.FontSize = 35; % Font size for tick labels
ax.LineWidth = 2;
axis equal
axis([-3 3 -3 3]) % Set axis limits to -3 to 3 for both x and y axes



figure()
scatter(points(:, 1), points(:, 2), [], abs(x_diff_1), 'filled');
colormap(gca, 'parula');
cb = colorbar;
set(cb, 'FontSize', 35);
cb.Position = [0.760, 0.1, 0.032 0.8]; % [left, bottom, width, height]
ax = gca; % Get current axes
ax.FontSize = 35; % Font size for tick labels
ax.LineWidth = 2;
axis equal
axis([-3 3 -3 3]) % Set axis limits to -3 to 3 for both x and y axes


% 
%%
%Using H_unif instead of H_opt 
%functions used are f(x, y) = 1 and f(x, y) = 4*sin(x)*sin(y)*e^(-(x^2+y^2))
N = size(points, 1);
D_x = N/(9*pi)*Q_x;
D_y = N/(9*pi)*Q_y;


%testing f(x, y) = 1
x_diff_1 = D_x*points(:, 1) - 1;
y_diff_1 = D_y*points(:, 2) - 1;



%testing f(x, y) = 4*sin(x)*sin(y)*e^(-(x^2+y^2))
x_diff_func = D_x*func(points(:, 1), points(:, 2)) - dx_func(points(:, 1), points(:, 2));
y_diff_func = D_y*func(points(:, 1), points(:, 2)) - dy_func(points(:, 1), points(:, 2));



%%
sqrt(sum((9*pi/N)*(x_diff_1.^2)))
sqrt(sum((9*pi/N)*(y_diff_1.^2)))
sqrt(sum((9*pi/N)*(x_diff_func.^2)))
sqrt(sum((9*pi/N)*(y_diff_func.^2)))



%%
function val = func(x, y)
    val = 4*sin(x).*sin(y).*exp(-(x.^2 + y.^2));
end

function val = dx_func(x, y)
    val = 4*sin(y).*exp(-(x.^2 + y.^2)).*(cos(x) - 2*x.*sin(x));
end

function val = dy_func(x, y)
    val = 4*sin(x).*exp(-(x.^2 + y.^2)).*(cos(y) - 2*y.*sin(y));
end