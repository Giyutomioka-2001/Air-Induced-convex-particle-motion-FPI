clc; clear; close all;

% Define grid sizes
N_values = [8, 16, 32, 64];

% Define functions to differentiate
functions = {@(x) 2*x, @(x) 5*x.^3};  
exact_derivatives = {@(x) 2*x, @(x) 15*x.^2};  
function_names = {'x^2', '5x^3'}; 

%% Loop over functions
for func_idx = 1:2
    f = functions{func_idx};
    f_prime_exact_func = exact_derivatives{func_idx};
    
    %% Figure 1: Padé vs. Exact
    figure;
    sgtitle(['Padé Scheme vs. Exact for ' function_names{func_idx}]);

    for idx = 1:length(N_values)
        N = N_values(idx);
        x = linspace(-1,1,N+1);
        h = x(2) - x(1);

        % Compute function values
        f_vals = f(x);

        % Compute derivatives
        f_prime_pade = pade_derivative(f_vals, h);

        % Compute analytical derivative
        f_prime_exact = f_prime_exact_func(x);

        % Subplot for each N value (Padé vs. Exact)
        subplot(2,2,idx);
        plot(x, f_prime_exact, 'k-', 'LineWidth', 1.5); hold on;
        plot(x, f_prime_pade, 'ro-', 'LineWidth', 1.2);
        title(['N = ' num2str(N)]);
        xlabel('x'); ylabel('f prime');
        legend('Exact', 'Padé');
        grid on;
    end

    %% Figure 2: Central Difference vs. Exact
    figure;
    sgtitle(['Central Difference vs. Exact for ' function_names{func_idx}]);

    for idx = 1:length(N_values)
        N = N_values(idx);
        x = linspace(-1,1,N+1);
        h = x(2) - x(1);

        % Compute function values
        f_vals = f(x);

        % Compute derivatives
        f_prime_cd = central_difference(f_vals, h);

        % Compute analytical derivative
        f_prime_exact = f_prime_exact_func(x);

        % Subplot for each N value (Central Difference vs. Exact)
        subplot(2,2,idx);
        plot(x, f_prime_exact, 'k-', 'LineWidth', 1.5); hold on;
        plot(x, f_prime_cd, 'b*-', 'LineWidth', 1.2);
        title(['N = ' num2str(N)]);
        xlabel('x'); ylabel('f prime');
        legend('Exact', 'Central Diff');
        grid on;
    end
end

hold off;

%% Function: Fourth-Order Padé Scheme using Tridiagonal System
function f_prime = pade_derivative(f, h)
    N = length(f);
    f_prime = zeros(size(f));
    A = zeros(N, N);
    b = zeros(N, 1);

    % Interior points
    for i = 2:N-1
        A(i, i-1) = 1;
        A(i, i) = 4;
        A(i, i+1) = 1;
        b(i) = (3 * (f(i+1) - f(i-1))) / h;
    end

    % Boundary conditions (Third-order accurate one-sided differences)
    A(1,1) = 1; A(1,2) = 2;
    b(1) = (-5/2*f(1) + 2*f(2) + f(3)/2) / h;
    A(N,N) = 1; A(N,N-1) = 2;
    b(N) = (5/2*f(N) - 2*f(N-1) - f(N-2)/2) / h;

    % Solve the tridiagonal system
    f_prime = A \ b;
end

%% Function: Second-Order Central Difference Scheme
function f_prime = central_difference(f, h)
    N = length(f);
    f_prime = zeros(size(f));

    % Interior points (Second-Order Central Difference)
    for i = 2:N-1
        f_prime(i) = (f(i+1) - f(i-1)) / (2*h);
    end

    % Boundary treatment (First-Order Forward and Backward Difference)
    f_prime(1) = (-3*f(1) + 4*f(2) - f(3)) / (2*h);
    f_prime(N) = (3*f(N) - 4*f(N-1) + f(N-2)) / (2*h);
end
