% MATLAB script to plot modified wavenumber k'h vs. kh for two numerical differentiation schemes

% Define kh range
kh = linspace(0, pi, 1000); % Generate values from 0 to pi

% Compute modified wavenumber for second-order central difference
kh_P_C = sin(kh);

% Compute modified wavenumber for fourth-order Pade scheme
kh_P_P = (6 * sin(kh)) ./ (4 + 2 * cos(kh));

% Exact wavenumber
kh_exact = kh;

% Plot the results
figure;
plot(kh, kh_exact, 'k--', 'LineWidth', 2); hold on;
plot(kh, kh_P_C, 'g', 'LineWidth', 2);
plot(kh, kh_P_P, 'black', 'LineWidth', 2);

% Labels and legend
xlabel('kh');
ylabel("k'h");
legend('Exact', '2nd-order Central Difference', '4th-order Pade Scheme');
title("Modified Wavenumber for exact vs numerical pade & central");
grid on;

% Display the plot
hold off;
