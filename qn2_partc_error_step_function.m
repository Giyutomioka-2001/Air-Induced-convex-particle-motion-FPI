clc; clear; close all;

% Define grid sizes (different N values)
N_values = [8, 16, 32, 64];
x_points = [-0.5, 0.4, 0.8];  % Selected interior points for error analysis
f = @(x) (x < 0);    % Step function
f_prime_exact = @(x) 0*x;  % Exact derivative is zero everywhere except at discontinuity

% Creating error array to store errors
errors_pade = zeros(length(N_values), length(x_points));
errors_cd = zeros(length(N_values), length(x_points));

% Loop through each grid size (N)
for idx = 1:length(N_values)
    N = N_values(idx);
    x = linspace(-1,1,N+1);   % Grid points in domain
    h = x(2) - x(1);         % Grid spacing
    
    % Compute function values
    f_vals = f(x);
    
    % Compute derivatives using Padé scheme (only for interior points)
    f_prime_pade = pade_derivative(f_vals, h);
    
    % Compute derivatives using Central Difference scheme
    f_prime_cd = central_difference(f_vals, h);
    
    % Compute exact derivatives at selected interior points
    for j = 1:length(x_points)
        point = x_points(j);
        
        % Ensure we pick only interior points for error calculation
        if point > x(1) && point < x(end)
            exact_value = f_prime_exact(point);
            
            % Find the index of the closest interior point in the grid
            [~, idx_point] = min(abs(x - point));
            
            % Compute errors at selected points
            errors_pade(idx, j) = abs(f_prime_pade(idx_point) - exact_value);
            errors_cd(idx, j) = abs(f_prime_cd(idx_point) - exact_value);
        else
            errors_pade(idx, j) = NaN;  % Ignore boundary points
            errors_cd(idx, j) = NaN;
        end
    end
end

% Plot the error as a function of number of intervals N
figure;
sgtitle('Error vs. Number of Intervals for Padé and Central Difference for Step Function');

% Plot Padé scheme error
subplot(2,1,1);
hold on;
for j = 1:length(x_points)
    if ~isnan(errors_pade(1, j))  % Plot only interior points
        loglog(N_values, errors_pade(:,j), 'o-', 'LineWidth', 1.5);
    end
end
xlabel('Number of intervals (N)');
ylabel('Absolute error');
legend(arrayfun(@(x) sprintf('Error at x = %.2f', x), x_points(~isnan(errors_pade(1,:))), 'UniformOutput', false), 'Location', 'NorthEast');
title('Padé Scheme Error');
grid on;

% Plot Central Difference scheme error
subplot(2,1,2);
hold on;
for j = 1:length(x_points)
    if ~isnan(errors_cd(1, j))  % Plot only interior points
        loglog(N_values, errors_cd(:,j), 'o-', 'LineWidth', 1.5);
    end
end
xlabel('Number of intervals (N)');
ylabel('Absolute error');
legend(arrayfun(@(x) sprintf('Error at x = %.2f', x), x_points(~isnan(errors_cd(1,:))), 'UniformOutput', false), 'Location', 'NorthEast');
title('Central Difference Scheme Error');
grid on;

%% Function: Fourth-Order Padé Scheme for Interior Points Only
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

    % Interior points only (excluding first and last points)
    for i = 2:N-1
        f_prime(i) = (f(i+1) - f(i-1)) / (2*h);
    end
end
