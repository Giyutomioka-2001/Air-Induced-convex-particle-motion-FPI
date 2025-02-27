clc; clear; close all;

% Define the domain for interpolation evaluation (fine grid)
x_domain = -1:0.01:1; % Fine grid with spacing 0.01

% Define the functions
sin_funcn = @(x) sin(x);% @ allows to calculate value of sin (x) at any given point without defining x and initialising it 
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
    sgtitle(['Lagrange Interpolation for ', funcn_title(f_index)]);

    % Loop over different numbers of sample points
    for index = 1:length(data_points)
        N = data_points(index); % Number of points

        % Sample points (equidistant)
        x_cal = linspace(-1, 1, N); % Equally spaced points
        y_cal = f(x_cal); % Compute function values at sample points

        % Perform Lagrange interpolation
        y_P = lagrange_interpolation(x_cal, y_cal, x_domain);

        % Plot the results
        subplot(2,2,index);% Multiple figures in single plot
        plot(x_domain, f(x_domain), 'green', 'LineWidth', 2.5); hold on;
        plot(x_domain, y_P, 'black', 'LineWidth', 1.5);
        plot(x_cal, y_cal, 'bo', 'MarkerFaceColor', 'black');
        title(sprintf('%d Sample Points', N));
        legend('Exact', 'Interpolation', 'Samples');
        grid on;
    end
end

% Lagrange Interpolation Function
function y_P = lagrange_interpolation(x_cal, y_cal, X)
    n = length(x_cal);
    y_P = zeros(size(X));% as summation of yi*Li so initialising it to 0 

    for i = 1:n
        L = ones(size(X));% as L will be multiplied so 1 not zero
        for j = 1:n
            if i ~= j
                L = L .* (X - x_cal(j)) / (x_cal(i) - x_cal(j));
            end
        end
        y_P = y_P + L * y_cal(i);
    end
end
