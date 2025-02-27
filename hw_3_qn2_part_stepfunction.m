clc; clear; close all;

% Define grid sizes
N_values = [8, 16, 32, 64];

for idx = 1:length(N_values)
    N = N_values(idx);
    fprintf('Tridiagonal Matrix for N = %d:\n', N);
    
    % Initialize tridiagonal matrix
    A = zeros(N, N);
    
    % Fill interior points
    for i = 2:N-1
        A(i, i-1) = 1;
        A(i, i) = 4;
        A(i, i+1) = 1;
    end
    
    % Boundary conditions
    A(1,1) = 1; A(1,2) = 2;
    A(N,N) = 1; A(N,N-1) = 2;
    
    % Display matrix
    disp(A);
    
    % Print matrix in formatted output
    fprintf('Tridiagonal Matrix for N = %d:\n', N);
    for i = 1:N
        for j = 1:N
            fprintf('%6.2f ', A(i,j));
        end
        fprintf('\n');
    end
    fprintf('\n');
end

% Define step function
step_function = @(x) (x < 0) .* 1 + (x > 0) .* 0;  % Step function

%% Plot Padé Scheme Results
figure;
sgtitle('Padé Scheme for Step Function');
for idx = 1:length(N_values)
    N = N_values(idx);
    x = linspace(-1,1,N+1);
    h = x(2) - x(1);

    % Compute function values
    f_vals = step_function(x);
    f_prime_pade = pade_derivative(f_vals, h);
    f_prime_exact = zeros(size(x));

    % Plot Padé Scheme
    subplot(2,2,idx);
    plot(x, f_prime_exact, 'k-', 'LineWidth', 1.5); hold on;
    plot(x, f_prime_pade, 'ro-', 'LineWidth', 1.2);
    title(['Padé Scheme: N = ' num2str(N)]);
    xlabel('x'); ylabel('f prime');
    legend('Exact', 'Padé');
    grid on;
end

%% Plot Central Difference Scheme Results
figure;
sgtitle('Central Difference Scheme for Step Function');
for idx = 1:length(N_values)
    N = N_values(idx);
    x = linspace(-1,1,N+1);
    h = x(2) - x(1);

    % Compute function values
    f_vals = step_function(x);
    f_prime_cd = central_difference(f_vals, h);
    f_prime_exact = zeros(size(x));

    % Plot Central Difference Scheme
    subplot(2,2,idx);
    plot(x, f_prime_exact, 'k-', 'LineWidth', 1.5); hold on;
    plot(x, f_prime_cd, 'b*-', 'LineWidth', 1.2);
    title(['Central Diff. Scheme: N = ' num2str(N)]);
    xlabel('x'); ylabel('f prime');
    legend('Exact', 'Central Diff.');
    grid on;
end

%% Fourth-Order Padé Scheme with Third-Order Boundary Closure
function f_prime = pade_derivative(f, h)
    N = length(f);
    f_prime = zeros(size(f));
    A = zeros(N, N);
    b = zeros(N, 1);

    % Interior points (Fourth-Order Padé)
    for i = 2:N-1
        A(i, i-1) = 1;
        A(i, i) = 4;
        A(i, i+1) = 1;
        b(i) = (3 * (f(i+1) - f(i-1))) / h;
    end

    % Boundary conditions (Third-order one-sided differences)
    A(1,1) = 1; A(1,2) = 2;
    b(1) = (-5/2*f(1) + 2*f(2) + f(3)/2) / h;
    A(N,N) = 1; A(N,N-1) = 2;
    b(N) = (5/2*f(N) - 2*f(N-1) - f(N-2)/2) / h;

    % Solve the system
    f_prime = A \ b;
end

%% Second-Order Central Difference Scheme
function f_prime = central_difference(f, h)
    N = length(f);
    f_prime = zeros(size(f));

    % Interior points
    for i = 2:N-1
        f_prime(i) = (f(i+1) - f(i-1)) / (2*h);
    end

    % Boundary conditions (First-order forward and backward difference)
    f_prime(1) = (-3*f(1) + 4*f(2) - f(3)) / (2*h);
    f_prime(N) = (3*f(N) - 4*f(N-1) + f(N-2)) / (2*h);
end
