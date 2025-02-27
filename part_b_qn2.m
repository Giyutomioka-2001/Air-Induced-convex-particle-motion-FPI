clc; clear; close all;

% Define the domain for interpolation evaluation (fine grid)
x_domain = -1:0.01:1; % Fine grid with spacing 0.01

% Define the functions
sin_funcn = @(x) sin(x);
step_funcn = @(x) (x < 0) * 1 + (x > 0) * 0;
runge_funcn = @(x) 1 ./ (1 + 25*x.^2);

% Store functions in a cell array for looping
func = {sin_funcn, step_funcn, runge_funcn};
funcn_title = ["sin(x)", "Step function", "Runge function"];

% Define number of sample points
data_points = [11, 21, 41, 81];

% Loop over each function
for f_index = 1:length(func)
    f = func{f_index}; % Select the function

    figure;
    sgtitle(['Spline Interpolation (Using Tridiagonal Matrix System) for ', funcn_title(f_index)]);

    % Loop for different number of sample points
    for index = 1:length(data_points)
        N = data_points(index); % Number of points

        % Sample points (equidistant)
        x_cal = linspace(-1, 1, N); % N points equally spaced between -1 to 1
        y_cal = f(x_cal); % Compute function values at sample points

        % Compute second derivatives using a tridiagonal system
        M = Second_prime(x_cal, y_cal);

        % Perform spline interpolation
        y_P = spline_interpolate(x_cal, y_cal, M, x_domain);

        % Plot the results
        subplot(2,2,index);
        plot(x_domain, f(x_domain), 'green', 'LineWidth', 1.5); hold on;
        plot(x_domain, y_P, 'black', 'LineWidth', 1.5);
        plot(x_cal, y_cal, 'bo', 'MarkerFaceColor', 'black');
        title(sprintf('%d Sample Points', N));
        legend('Exact', 'Interpolation', 'Samples');
        grid on;
    end
end

% Function to compute second derivatives using a tridiagonal system
function M = Second_prime(x_cal, y_cal)
    n = length(x_cal);
    h = diff(x_cal);
    A = zeros(n, n);
    rhs = zeros(n, 1);
    
    for i = 2:n-1
        A(i, i-1) = h(i-1);
        A(i, i) = 2 * (h(i-1) + h(i));
        A(i, i+1) = h(i);
        rhs(i) = (6/h(i)) * (y_cal(i+1) - y_cal(i)) - (6/h(i-1)) * (y_cal(i) - y_cal(i-1));
    end
    
    A(1,1) = 1; % Free run-out condition
    A(n,n) = 1; % Free run-out condition
    
    % Solve the system
    M = A \ rhs;
end

% Function to evaluate spline interpolation
function y_interp = spline_interpolate(x_cal, y_cal, M, x_interp)
    n = length(x_cal) - 1;
    h = diff(x_cal);
    y_interp = zeros(size(x_interp));
    
    for i = 1:n
        idx = (x_interp >= x_cal(i)) & (x_interp <= x_cal(i+1));
        xi = x_cal(i);
        xi1 = x_cal(i+1);
        yi = y_cal(i);
        yi1 = y_cal(i+1);
        hi = xi1 - xi;
        
        A = (xi1 - x_interp(idx)) / hi;
        B = (x_interp(idx) - xi) / hi;
        C = (A.^3 - A) * (hi^2) / 6;
        D = (B.^3 - B) * (hi^2) / 6;
        
        y_interp(idx) = A .* yi + B .* yi1 + C .* M(i) + D .* M(i+1);
    end
end
