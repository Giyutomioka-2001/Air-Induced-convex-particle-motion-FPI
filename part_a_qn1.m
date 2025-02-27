clc; clear; close all;

% Define the domain
x_domain = 0:0.01:1; % Fine grid for plotting
sample_pts = [4, 8, 16]; % Different numbers of sample points

% Define polynomial degrees:
poly_degrees = { [2, 4, 6], [4, 8, 16], [7, 15, 19] }; % Less, Equal, Greater

figure;
for index = 1:length(sample_pts)
    num_points = sample_pts(index); % Number of sampling points
    degrees = poly_degrees{index}; % Get corresponding polynomial degrees

    for d_index = 1:length(degrees)
        n = degrees(d_index); % Polynomial degree

        % Define function to interpolate
        f = @(x) x.^n; % Function f(x) = x^n

        x_cal = linspace(0, 1, num_points); % Equidistant interpolation points
        y_cal = f(x_cal); % Compute function values at sample points

        % Perform Lagrange interpolation
        y_P = lagrange_interpolation(x_cal, y_cal, x_domain);

        % Plot results
        subplot(length(sample_pts), 3, (index-1)*3 + d_index);
        plot(x_domain, f(x_domain), 'green', 'LineWidth', 2.5); hold on;
        plot(x_domain, y_P, 'black', 'LineWidth', 1.5);
        plot(x_cal, y_cal, 'bo', 'MarkerFaceColor', 'black','MarkerSize', 4);
        title(sprintf('$x^{%d}$ with %d points', n, num_points), 'Interpreter', 'latex');
        legend('Exact', 'Interpolated', 'Sample data points');
        grid on;
    end
end

sgtitle('Lagrange Interpolation for f(x) = x^n');

% Lagrange Interpolation Function
function y_P = lagrange_interpolation(x_cal, y_cal, X)
    n = length(x_cal);
    y_P = zeros(size(X));

    for i = 1:n
        L = ones(size(X));
        for j = 1:n
            if i ~= j
                L = L .* (X - x_cal(j)) / (x_cal(i) - x_cal(j));
            end
        end
        y_P = y_P + L * y_cal(i);
    end
end
