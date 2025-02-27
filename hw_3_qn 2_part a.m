clc; clear; close all;

% Define grid sizes
N_values = [8, 16, 32, 64];

% Define functions to differentiate
f1 = @(x) sin(pi*x);  % Function 1: sin(pi*x)
f2 = @(x) (x < 0) .* 1; % Function 2: Step function

% Loop over different grid sizes
for N = N_values
    % Discretize domain
    x = linspace(-1,1,N+1);
    h = x(2) - x(1);
    
    % Compute function values
    f1_vals = f1(x);
    f2_vals = f2(x);

    % Compute first derivatives
    dfdx_pade_f1 = pade_derivative(f1_vals, h);
    dfdx_cd_f1 = central_difference(f1_vals, h);

    dfdx_pade_f2 = pade_derivative(f2_vals, h);
    dfdx_cd_f2 = central_difference(f2_vals, h);

    % Compute analytical derivatives
    dfdx_exact_f1 = pi * cos(pi * x);
    dfdx_exact_f2 = zeros(size(x)); % Discontinuous

    % Plot results for f1 = sin(pi*x)
    figure;
    plot(x, dfdx_exact_f1, 'k-', 'LineWidth', 1.5); hold on;
    plot(x, dfdx_pade_f1, 'ro-', 'LineWidth', 1.2);
    plot(x, dfdx_cd_f1, 'b*-', 'LineWidth', 1.2);
    xlabel('x'); ylabel('df/dx');
    title(['Derivative of sin(\pi x) for N = ' num2str(N)]);
    legend('Exact', 'Padé', 'Central Diff');
    grid on;

    % Plot error for f1
    figure;
    plot(x, abs(dfdx_pade_f1 - dfdx_exact_f1), 'ro-', 'LineWidth', 1.2); hold on;
    plot(x, abs(dfdx_cd_f1 - dfdx_exact_f1), 'b*-', 'LineWidth', 1.2);
    xlabel('x'); ylabel('Error');
    title(['Error in Derivative of sin(\pi x) for N = ' num2str(N)]);
    legend('Padé Error', 'Central Diff Error');
    grid on;
end

% Function: Fourth-Order Padé Scheme
function dfdx = pade_derivative(f, h)
    N = length(f);
    dfdx = zeros(size(f));

    % Interior points (Padé scheme)
    for i = 2:N-1
        dfdx(i) = (3/2) * (f(i+1) - f(i-1)) / (2*h);
    end

    % Boundary treatment (Third-order accurate one-sided differences)
    dfdx(1) = (-3*f(1) + 4*f(2) - f(3)) / (2*h);
    dfdx(N) = (3*f(N) - 4*f(N-1) + f(N-2)) / (2*h);
end

% Function: Second-Order Central Difference Scheme
function dfdx = central_difference(f, h)
    N = length(f);
    dfdx = zeros(size(f));

    % Interior points (Central Difference)
    for i = 2:N-1
        dfdx(i) = (f(i+1) - f(i-1)) / (2*h);
    end

    % Boundary treatment (First-order forward and backward difference)
    dfdx(1) = (-3*f(1) + 4*f(2) - f(3)) / (2*h);
    dfdx(N) = (3*f(N) - 4*f(N-1) + f(N-2)) / (2*h);
end
