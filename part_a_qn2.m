clc; clear; close all;

%% Part (a) - Spline Interpolation for Polynomial Test Functions (x^n)
% Define the domain
x_domain = linspace(0, 1, 500); % Fine grid for plotting
sample_pts = [4, 8, 16]; % Different numbers of sample points

% Define polynomial degrees (Less than, Equal, Greater than sample points)
poly_degrees = { [2, 4, 6], [6, 8, 10], [14, 16, 18] }; 

figure;
for index = 1:length(sample_pts)
    num_points = sample_pts(index); % Get number of sample points
    degrees = poly_degrees{index}; % Get corresponding polynomial degrees

    for d_index = 1:length(degrees)
        n = degrees(d_index); % Polynomial degree

        % Define function
        f = @(x) x.^n;

        % Sample points
        x_cal = linspace(0, 1, num_points);
        y_cal = f(x_cal);

        % Perform Spline Interpolation using MATLAB's built-in function
        y_S = spline(x_cal, y_cal, x_domain);

        % Plot results
        subplot(length(sample_pts), 3, (index-1)*3 + d_index);
        plot(x_domain, f(x_domain), 'green', 'LineWidth', 2.5); hold on;
        plot(x_domain, y_S, 'black', 'LineWidth', 1.5);
        plot(x_cal, y_cal, 'bo', 'MarkerFaceColor', 'black', 'MarkerSize', 4);
        title(sprintf('$x^{%d}$ with %d points', n, num_points), 'Interpreter', 'latex');
        legend('Exact', 'Spline Interpolation', 'Samples', 'FontSize', 8);
        grid on;
    end
end

sgtitle('Spline Interpolation for f(x) = x^n using tridiagonal system of matrix');

%% Part (b) - Spline Interpolation for sin(x), Step Function, and Runge Function
x_domain = -1:0.01:1; % Fine grid for evaluation

% Define functions
sin_funcn = @(x) sin(x);
step_funcn = @(x) (x < 0) * 1 + (x > 0) * 0;
runge_funcn = @(x) 1 ./ (1 + 25*x.^2);

func = {sin_funcn, step_funcn, runge_funcn};
funcn_title = ["sin(x)", "Step function", "Runge function"];

% Define number of sample points
data_points = [11, 21, 41, 81];

% Loop over each function
for f_index = 1:length(func)
    f = func{f_index};

    figure;
    sgtitle(['Spline Interpolation for ', funcn_title(f_index)]);

    for index = 1:length(data_points)
        N = data_points(index);

        % Sample points
        x_cal = linspace(-1, 1, N);
        y_cal = f(x_cal);

        % Perform Spline Interpolation
        y_S = spline(x_cal, y_cal, x_domain);

        % Plot results
        subplot(2,2,index);
        plot(x_domain, f(x_domain), 'green', 'LineWidth', 2.5); hold on;
        plot(x_domain, y_S, 'black', 'LineWidth', 1.5);
        plot(x_cal, y_cal, 'bo', 'MarkerFaceColor', 'black', 'MarkerSize', 4);
        title(sprintf('%d Sample Points', N));
        legend('Exact', 'Spline Interpolation', 'Samples', 'FontSize', 8);
        grid on;
    end
end

