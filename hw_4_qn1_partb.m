clc; clear; close all;

a = 0.98;   
N_nonuniform = 32;  % Non-uniform grid points

% === Non-uniform grid ===
j = 0:N_nonuniform;
xi = -1 + (2 * j) / N_nonuniform; 
x_nonuniform = (1/a) * tanh(xi * atanh(a)); % Non-uniform grid points

% Function and exact derivative
f = 1 - x_nonuniform.^8;
df_exact = -8 * x_nonuniform.^7;

% Numerical derivative for non-uniform grid (central difference)
df_num_nonuniform = zeros(1, N_nonuniform+1);

for j = 2:N_nonuniform
    h_j = x_nonuniform(j) - x_nonuniform(j-1);
    h_j1 = x_nonuniform(j+1) - x_nonuniform(j);
    df_num_nonuniform(j) = (f(j+1) - f(j-1)) / (h_j1 + h_j);
end

df_num_nonuniform(1) = (f(2) - f(1)) / (x_nonuniform(2) - x_nonuniform(1));
df_num_nonuniform(N_nonuniform+1) = (f(N_nonuniform+1) - f(N_nonuniform)) / (x_nonuniform(N_nonuniform+1) - x_nonuniform(N_nonuniform));

% Compute maximum error for non-uniform grid
max_error_nonuniform = max(abs(df_num_nonuniform - df_exact));

% === Uniform grid comparison ===
N_uniform = 32; % Start with same number of points
converged = false;

while ~converged
    x_uniform = linspace(-1, 1, N_uniform+1);  % Uniform grid
    f_uniform = 1 - x_uniform.^8;
    df_exact_uniform = -8 * x_uniform.^7;
    
    df_num_uniform = zeros(1, N_uniform+1);
    
    % Central difference for uniform grid
    for j = 2:N_uniform
        h = x_uniform(j+1) - x_uniform(j);
        df_num_uniform(j) = (f_uniform(j+1) - f_uniform(j-1)) / (2 * h);
    end
    
    df_num_uniform(1) = (f_uniform(2) - f_uniform(1)) / (x_uniform(2) - x_uniform(1));
    df_num_uniform(N_uniform+1) = (f_uniform(N_uniform+1) - f_uniform(N_uniform)) / (x_uniform(N_uniform+1) - x_uniform(N_uniform));
    
    % Compute max error for uniform grid
    max_error_uniform = max(abs(df_num_uniform - df_exact_uniform));
    
    % Check if uniform grid error is ≤ non-uniform grid error
    if max_error_uniform <= max_error_nonuniform
        converged = true;
    else
        N_uniform = N_uniform + 2; % Increase resolution
    end
end

% Display results
fprintf('Max error (non-uniform grid, N = %d): %e\n', N_nonuniform, max_error_nonuniform);
fprintf('Required N for uniform grid: %d (Max error = %e)\n', N_uniform, max_error_uniform);
