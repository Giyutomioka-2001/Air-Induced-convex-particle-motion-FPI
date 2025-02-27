clc; clear; close all;

% Define grid sizes
N_values = [8, 16, 32, 64];

% Define polynomial test functions
test_functions = {
    @(x) ones(size(x)), @(x) x, @(x) x.^2, @(x) x.^3
};
test_derivatives = {
    @(x) zeros(size(x)), @(x) ones(size(x)), @(x) 2*x, @(x) 3*x.^2
};

% Run tests
fprintf('Testing Padé scheme on polynomial functions...\n');

for idx = 1:length(N_values)
    N = N_values(idx);
    x = linspace(-1,1,N+1);
    h = x(2) - x(1);

    for j = 1:length(test_functions)
        f = test_functions{j};
        df_exact = test_derivatives{j};

        % Compute function values
        f_vals = f(x);

        % Compute Padé derivative
        f_prime_pade = pade_derivative(f_vals, h);

        % Compute exact derivative
        f_prime_exact = df_exact(x);

        % Check error (should be close to zero for polynomials of order 3 or less)
        error = max(abs(f_prime_pade - f_prime_exact));
        
        fprintf('N = %d | Polynomial Order = %d | Max Error = %.2e\n', N, j-1, error);
        
        % Assert if error is above machine precision
        assert(error < 1e-12, 'Test failed for polynomial of order %d', j-1);
    end
end

fprintf('All polynomial tests passed! Padé scheme is correctly implemented.\n');
