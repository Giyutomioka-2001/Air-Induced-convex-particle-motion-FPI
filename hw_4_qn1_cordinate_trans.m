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
f = 1 - x.^8;  % f(x) = 1 - x^8
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
xlabel('x'); ylabel("f'(x)");
legend('dF exact', 'dF- Central Difference');
title('F prime using Central Difference for Non-Uniform Grid');
grid on;