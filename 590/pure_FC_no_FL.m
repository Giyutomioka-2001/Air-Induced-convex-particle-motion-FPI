clc; clear; close all;

%% Simulation Parameters
N_particles = 720;                 % Number of pentagon particles
pent_edge = 0.005;                % Edge length in meters (5 mm)
bounds = [0, 0.32; 0, 0.32];      % 32 cm x 32 cm domain
m = 1;                            % Uniform mass
cr = 1;                           % Coefficient of restitution
timeStep = 0.1;                  % Time step (can be tuned)
numSteps = 1000;                   % Total steps

% Particle area (regular pentagon)
A_pentagon = (5/4) * pent_edge^2 / tan(pi/5);  % Formula for regular pentagon area

%% Create Cartesian Grid
grid_dx = 0.02;  % 2 cm
grid_dy = 0.02;
x_grid = bounds(1,1):grid_dx:bounds(1,2);
y_grid = bounds(2,1):grid_dy:bounds(2,2);
nx = length(x_grid) - 1;
ny = length(y_grid) - 1;

%% Initialize Particles
particles = struct();
positions = [];

for i = 1:N_particles
    while true
        pos = (bounds(:,2) - bounds(:,1)) .* rand(2,1) + bounds(:,1);
        if isempty(positions) || all(vecnorm(positions - pos, 2, 1) > 2 * pent_edge)
            positions = [positions, pos];
            break;
        end
    end
    
    % Pentagon vertices centered at origin
    angle = linspace(0, 2*pi, 6);
    x = pent_edge * cos(angle);
    y = pent_edge * sin(angle);
    
    particles(i).vertices = [x; y];
    particles(i).position = pos;
    particles(i).velocity = 0.1 * (2 * rand(2,1) - 1);  % Small random velocity
end

%% Simulation Loop
figure;
for step = 1:numSteps
    clf; hold on; axis equal;
    xlim(bounds(1, :)); ylim(bounds(2, :));
    title(['Step: ', num2str(step)]);

    % Initialize area fraction grid
    area_fraction = zeros(nx, ny);

    % Update positions and compute area fractions
    for i = 1:N_particles
        % Update position
        particles(i).position = particles(i).position + particles(i).velocity * timeStep;

        % Boundary conditions
        for dim = 1:2
            if particles(i).position(dim) < bounds(dim, 1) || particles(i).position(dim) > bounds(dim, 2)
                particles(i).position(dim) = max(min(particles(i).position(dim), bounds(dim, 2)), bounds(dim, 1));
                particles(i).velocity(dim) = -cr * particles(i).velocity(dim);
            end
        end

        % Draw particle
        verticesWorld = particles(i).vertices + particles(i).position;
        fill(verticesWorld(1, :), verticesWorld(2, :), 'b', 'FaceAlpha', 0.5);

        % Determine grid cell the particle center is in
        x_idx = find(x_grid <= particles(i).position(1), 1, 'last');
        y_idx = find(y_grid <= particles(i).position(2), 1, 'last');

        if x_idx < nx+1 && y_idx < ny+1
            A_cell = grid_dx * grid_dy;
            area_fraction(x_idx, y_idx) = area_fraction(x_idx, y_idx) + A_pentagon / A_cell;
        end
    end

    % Visualize grid
    for x = x_grid
        plot([x x], bounds(2,:), 'k:', 'LineWidth', 0.5);
    end
    for y = y_grid
        plot(bounds(1,:), [y y], 'k:', 'LineWidth', 0.5);
    end

    % Display area fraction as red text
    for ix = 1:nx
        for iy = 1:ny
            cx = (x_grid(ix) + x_grid(ix+1)) / 2;
            cy = (y_grid(iy) + y_grid(iy+1)) / 2;
            val = area_fraction(ix, iy);
            if val > 0
                text(cx, cy, sprintf('%.2f', val), 'Color', 'r', 'FontSize', 8, ...
                    'HorizontalAlignment', 'center');
            end
        end
    end

    % Collision handling
    for i = 1:N_particles
        for j = i+1:N_particles
            dist = norm(particles(i).position - particles(j).position);
            min_dist = 2.0 * pent_edge;

            if dist < min_dist
                normal = (particles(i).position - particles(j).position) / (dist + eps);
                correction = (min_dist - dist) / 2;
                particles(i).position = particles(i).position + correction * normal;
                particles(j).position = particles(j).position - correction * normal;

                v_rel = dot(particles(i).velocity - particles(j).velocity, normal);
                if v_rel < 0
                    J_c = -(1 + cr) * v_rel / (1/m + 1/m);
                    particles(i).velocity = particles(i).velocity + (J_c / m) * normal;
                    particles(j).velocity = particles(j).velocity - (J_c / m) * normal;
                end
            end
        end
    end

    drawnow;
    pause(0.01);
end
