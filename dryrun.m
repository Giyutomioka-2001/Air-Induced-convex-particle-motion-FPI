% Define the number of intervals
intervals = [8, 16, 32];

% Define selected points where error will be measured
selected_points = [-0.5, 0, 0.5]; 

% Initialize error matrix
errors = zeros(length(selected_points), length(intervals));

% Loop over selected points
for i = 1:length(selected_points)
    
    % Loop over different grid resolutions
    for j = 1:length(intervals)
        N = intervals(j);  % Number of intervals
        x = linspace(-1, 1, N+1);  % Grid points

        % Define the test function and its exact derivative
        f_test = x.^4;
        f_test_exact = 4 * x.^3;

        % Compute the numerical derivative using the Pade scheme
        f_test_pade = pade_derivative(x, f_test);

        % Find the closest grid point to the selected point
        [~, idx] = min(abs(x - selected_points(i)));

        % Compute the error at the selected point
        errors(i, j) = abs(f_test_pade(idx) - f_test_exact(idx));
    end
end

% Display errors
disp('Errors at selected points for different grid resolutions:')
disp(errors)
