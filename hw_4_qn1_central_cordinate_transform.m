clc; clear; close all;
a = 0.98;   
N = 32;    

% non-uniform grid
j = 0:N;
xi = -1 + (2 * j) / N; 
x = (1/a) * tanh(xi * atanh(a)); % non-uniform grid points

% Function 
f = 1 - x.^8;  % f(x) = 1 - x^8
df_exact = -8 * x.^7;  % Analytical derivative f'(x)

% Computing numerical derivative using central difference
df_num = zeros(1, N+1); % Initializing derivative array

for j = 2:N
    h_j   = x(j) - x(j-1);   % h_j = x_j - x_(j-1)
    h_j1  = x(j+1) - x(j);   % h_(j+1) = x_(j+1) - x_j
    df_num(j) = (f(j+1) - f(j-1)) / (h_j1 + h_j); % using central difference formula
end

% forward and backward difference at boundaries
df_num(1) = (f(2) - f(1)) / (x(2) - x(1)); %fwd
df_num(N+1) = (f(N+1) - f(N)) / (x(N+1) - x(N));  %bwd

% Plots
figure;
plot(x, df_exact, 'k-', 'LineWidth', 2); hold on;
plot(x, df_num, 'go--', 'MarkerFaceColor', 'b');
xlabel('x'); ylabel("f'(x)");
legend('dF exact', 'dF- Central Difference');
title('F prime using Central Difference for Non-Uniform Grid');
grid on;
%
figure;
plot(x, df_exact, 'k-', 'LineWidth', 2); hold on;
plot(x, df_num, 'go--', 'MarkerFaceColor', 'b');

% to visualize non-uniform spacing
y_marker = min(df_exact) - (max(df_exact) - min(df_exact)); % Adjust marker level
scatter(x, y_marker * ones(size(x)), 15, 'r', 'filled'); % Small red dots on x-axis

xlabel('x'); ylabel("f'(x)");
legend('dF exact', 'dF- Central Difference');
title('F prime using Central Difference for Non-Uniform Grid');
grid on;


%% Central difference using cordinate transformation
clc; clear; close all;

% Given parameters
a = 0.98;   % Grid parameter
N = 32;     % Number of grid points

% Define ξ_j values
j = 0:N;
xi = -1 + (2 * j) / N;  % ξ_j = -1 + 2j/N

% Compute the non-uniform grid x_j
x = (1/a) * tanh(xi * atanh(a));

% Define function f(x) and its exact derivative
f = 1-x.^8;  % f(x) = 1 - x^8
df_exact = -8 * x.^7;  % Analytical derivative f'(x)

% Compute df/dξ using central difference
df_dxi = zeros(1, N+1);  % Initialize df/dξ array

for j = 2:N
    dxi_j = xi(j) - xi(j-1);    % Δξ_j = ξ_j - ξ_(j-1)
    dxi_j1 = xi(j+1) - xi(j);   % Δξ_(j+1) = ξ_(j+1) - ξ_j
    df_dxi(j) = (f(j+1) - f(j-1)) / (dxi_j1 + dxi_j);
end

% Use forward and backward difference at boundaries
df_dxi(1) = (f(2) - f(1)) / (xi(2) - xi(1));
df_dxi(N+1) = (f(N+1) - f(N)) / (xi(N+1) - xi(N));

% Compute dξ/dx
dxi_dx = (a / atanh(a)) * cosh(xi * atanh(a)).^2;

% Compute final derivative f'(x)
df_num = df_dxi .* dxi_dx;

% Plot results

figure;
plot(x, df_exact, 'k-', 'LineWidth', 2); hold on;
plot(x, df_num, 'go--', 'MarkerFaceColor', 'b');

% to visualize non-uniform spacing
y_marker = min(df_exact) - (max(df_exact) - min(df_exact)); % Adjust marker level
scatter(x, y_marker * ones(size(x)), 15, 'black', 'filled'); % Small red dots on x-axis

xlabel('x'); ylabel("f'(x)");
legend('dF exact', 'dF- Central Difference');
title('F prime using Central Difference for Non-Uniform Grid');
grid on;
