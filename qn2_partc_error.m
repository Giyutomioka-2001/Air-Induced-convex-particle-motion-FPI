clc; clear; close all;

% Define grid sizes (different N values)
N_values = [8, 16, 32, 64];
x_points = [-0.5, 0.4, 0.8];  % Interior points for error analysis
f = @(x) sin(pi*x);    % Function: sin(pi*x)
f_prime_exact = @(x) pi * cos(pi*x);  % Exact derivative: pi * cos(pi*x)

% Error storage
errors_pade = zeros(length(N_values), length(x_points));
errors_cd = zeros(length(N_values), length(x_points));

% Loop through each grid size (N)
for idx = 1:length(N_values)
    N = N_values(idx);
    x = linspace(-1,1,N+1);   % Grid points in domain
    h = x(2) - x(1);         % Grid spacing
    
    % Compute function values
    f_vals = f(x);
    
    % Compute derivatives using Padé scheme (interior points only)
    f_prime_pade = pade_derivative(f_vals, h);
    
    % Compute derivatives using Central Difference scheme
    f_prime_cd = central_difference(f_vals, h);
    
    % Compute exact derivatives at selected interior points
    for j = 1:length(x_points)
        point = x_points(j);
        
        % Find the closest grid index for the selected interior points
        [~, idx_point] = min(abs(x - point));
        
        % Ensure the selected point is an actual interior point
        if idx_point > 1 && idx_point < N+1
            exact_value = f_prime_exact(x(idx_point));
            errors_pade(idx, j) = abs(f_prime_pade(idx_point) - exact_value);
            errors_cd(idx, j) = abs(f_prime_cd(idx_point) - exact_value);
        else
            errors_pade(idx, j) = NaN;
            errors_cd(idx, j) = NaN;
        end
    end
end

% Plot Error
figure;
sgtitle('Error vs. Number of Intervals for Padé and Central Difference');

% Padé scheme error plot
subplot(2,1,1);
hold on;
for j = 1:length(x_points)
    loglog(N_values, errors_pade(:,j), 'o-', 'LineWidth', 1.5);
end
xlabel('Number of intervals (N)');
ylabel('Absolute error');
legend(arrayfun(@(x) sprintf('x = %.2f', x), x_points, 'UniformOutput', false), 'Location', 'NorthEast');
title('Padé Scheme Error');
grid on;

% Central Difference error plot
subplot(2,1,2);
hold on;
for j = 1:length(x_points)
    loglog(N_values, errors_cd(:,j), 'o-', 'LineWidth', 1.5);
end
xlabel('Number of intervals (N)');
ylabel('Absolute error');
legend(arrayfun(@(x) sprintf('x = %.2f', x), x_points, 'UniformOutput', false), 'Location', 'NorthEast');
title('Central Difference Scheme Error');
grid on;

%% Fourth-Order Padé Scheme Function
function f_prime = pade_derivative(f, h)
    N = length(f);
    f_prime = zeros(N,1);
    A = diag(4*ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
    b = zeros(N,1);
    
    for i = 2:N-1
        b(i) = (3 * (f(i+1) - f(i-1))) / h;
    end
    
    % Boundary Conditions (One-Sided Differences)
    b(1) = (-5/2*f(1) + 2*f(2) + f(3)/2) / h;
    b(N) = (5/2*f(N) - 2*f(N-1) - f(N-2)/2) / h;
    
    f_prime = A \ b; % Solve system
end

%% Second-Order Central Difference Function
function f_prime = central_difference(f, h)
    N = length(f);
    f_prime = zeros(N,1);
    
    for i = 2:N-1
        f_prime(i) = (f(i+1) - f(i-1)) / (2*h);
    end
end
