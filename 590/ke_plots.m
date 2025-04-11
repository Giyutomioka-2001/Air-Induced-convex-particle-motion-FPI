%% KE Elastic vs Inelastic
clc; clear; close all;

%Simulation Parameters
N_particles = 100;
pent_edge = 0.005;              % Edge length in meters (5 mm)
bounds = [0, 0.32; 0, 0.32];    % 32 cm x 32 cm domain
m = 1;                          % Uniform mass
timeStep = 0.01;                % Time step
numSteps = 500;

% Particle area
A_pentagon = (5/4) * pent_edge^2 / tan(pi/5);

% Grid
grid_dx = 0.02;
grid_dy = 0.02;
x_grid = bounds(1,1):grid_dx:bounds(1,2);
y_grid = bounds(2,1):grid_dy:bounds(2,2);
nx = length(x_grid) - 1;
ny = length(y_grid) - 1;

% Run for both cr values
cr_list = [1, 0.8];
KE_data = zeros(numSteps, length(cr_list));

for mode = 1:length(cr_list)
    cr = cr_list(mode);

    % Initialize Particles
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

        angle = linspace(0, 2*pi, 6);
        x = pent_edge * cos(angle);
        y = pent_edge * sin(angle);
        particles(i).vertices = [x; y];
        particles(i).position = pos;
        particles(i).velocity = 0.1 * (2 * rand(2,1) - 1);
    end

    % Simulation Loop
    for step = 1:numSteps
        area_fraction = zeros(nx, ny);
        KE = 0; % Kinetic energy per timestep

        for i = 1:N_particles
            particles(i).position = particles(i).position + particles(i).velocity * timeStep;

            for dim = 1:2
                if particles(i).position(dim) < bounds(dim, 1) || particles(i).position(dim) > bounds(dim, 2)
                    particles(i).position(dim) = max(min(particles(i).position(dim), bounds(dim, 2)), bounds(dim, 1));
                    particles(i).velocity(dim) = -cr * particles(i).velocity(dim);
                end
            end

            % KE tracking
            KE = KE + 0.5 * m * norm(particles(i).velocity)^2;

            % Area fraction grid update
            x_idx = find(x_grid <= particles(i).position(1), 1, 'last');
            y_idx = find(y_grid <= particles(i).position(2), 1, 'last');
            if x_idx < nx+1 && y_idx < ny+1
                A_cell = grid_dx * grid_dy;
                area_fraction(x_idx, y_idx) = area_fraction(x_idx, y_idx) + A_pentagon / A_cell;
            end
        end

        % Store KE
        KE_data(step, mode) = KE;

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
    end
end

% Plot KE vs Time
figure;
plot(1:numSteps, KE_data(:,1), 'b-', 'LineWidth', 2); hold on;
plot(1:numSteps, KE_data(:,2), 'r-', 'LineWidth', 2);
xlabel('Time Step'); ylabel('Total Kinetic Energy');
legend('Elastic (cr = 1)', 'Inelastic (cr = 0.8)', 'Location', 'northeast');
title('Kinetic Energy vs Time for Elastic and Inelastic Collisions');
grid on;
