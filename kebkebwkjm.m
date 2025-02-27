% Define the range for x
x = linspace(0, 2*pi, 100);

% Compute y values for sin(x) and cos(x)
y1 = sin(x);
y2 = cos(x);

% Create the plot
figure;
plot(x, y1, 'g', 'LineWidth', 2); % Plot sin(x) in red
hold on;
plot(x, y2, 'c', 'LineWidth', 2); % Plot cos(x) in blue
hold off;

% Add labels, title, and legend
xlabel('x (radians)');
ylabel('y');
title('Graphs of y = sin(x) and y = cos(x)');
legend('y = sin(x)', 'y = cos(x)', 'Location', 'best');

% Customize the grid and axes
grid on;
axis([0, 2*pi, -1.5, 1.5])
