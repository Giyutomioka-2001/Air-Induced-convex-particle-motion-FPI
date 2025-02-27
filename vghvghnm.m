% Define the range for t
t = linspace(0, 4*pi, 100);

% Parametric functions for the helix
x = cos(t);
y = sin(t);
z = t;

% Create the 3D plot
figure;
plot3(x, y, z, 'r-', 'LineWidth', 2); % Plot the helix with a red solid line
grid on;

% Add labels and title
xlabel('x = cos(t)');
ylabel('y = sin(t)');
zlabel('z = t');
title('3D Helix');

% Customize the view
view(45, 30); % Set the view angle for better visualization
axis tight;   % Fit the axes tightly around the curve
