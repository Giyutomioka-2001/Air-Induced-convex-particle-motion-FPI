% MATLAB Code for Simulating Pentagon-Shaped Particles in a Circular Ring (2D)
% Particles in circular motion
clc; clear; close all;

%% Simulation Parameters
N_particles = 15; % Number of pentagon particles
radius = 20; % Radius of circular ring in cm
pent_edge = 1.0 ; % Edge length of the pentagon in cm
bounds = [-10, 10; -10, 10]; % Simulation box bounds [x_min, x_max; y_min, y_max]
m = 1; % Uniform m for all particles
cr = 0.9; % Higher coefficient of restitution for more collisions

Cl = 0.8; % Assumed moderate lift coefficient
V_air = 1.0; % Moderate air velocity
rho_air = 1.2; % Air density SI units

% Compute apothem
apothem = pent_edge / (2 * tan(pi/5)); % Apothem for a regular pentagon
A_pentagon = (5/2) * pent_edge * apothem; % Approx. area of pentagon

timeStep = 0.2; % Time step for simulation
numSteps = 2500; % Number of simulation steps

%% Define Grid Parameters
Nx = 20; % Number of grid cells in x-direction
Ny = 20; % Number of grid cells in y-direction
Lx = bounds(1, 2) - bounds(1, 1); % Length of simulation domain in x-direction
Ly = bounds(2, 2) - bounds(2, 1); % Length of simulation domain in y-direction
dx = Lx / Nx; % Grid cell size in x-direction
dy = Ly / Ny; % Grid cell size in y-direction

% Initialize grid for area fraction
area_fraction_grid = zeros(Nx, Ny);

%% Compute Area Fraction (φ_A)
phi_A = (N_particles * (5/2)* pent_edge * apothem) / (pi * radius^2); % Area fraction calculation
fprintf('Area fraction (φ_A): %.4f\n', phi_A); % Print area fraction

%% Initialize Particles
particles = struct();
positions = [];
for i = 1:N_particles
    while true
        pos = 8 * rand(2,1) - [1;5]; % Random position
        if isempty(positions) || all(vecnorm(positions - pos, 2, 1) > 1.0 * pent_edge)
            positions = [positions, pos];
            break;
        end
    end
    
    % Define pentagon vertices (regular pentagon)
    angle = linspace(0, 2*pi, 6); % One extra point to close the pentagon
    x = pent_edge * cos(angle);
    y = pent_edge * sin(angle);
    particles(i).vertices = [x; y];
    particles(i).position = pos;
    particles(i).velocity = rand(2,1) * 2 - 1; % Random initial velocity
    particles(i).m = m;
    particles(i).lift_force = [0; 0]; % Initialize lift force
end

%% Simulation Loop
figure;
for step = 1:numSteps
    clf;
    hold on;
    axis equal;
    xlim(bounds(1, :));
    ylim(bounds(2, :));
    
    % Reset area fraction grid
    area_fraction_grid = zeros(Nx, Ny);
    
    for i = 1:N_particles
        % Compute lift force in the x-y plane
        v_mag = norm(particles(i).velocity);
        lift_magnitude = Cl * (rho_air * V_air^2 / 2) * A_pentagon;
        
        % Lift direction is perpendicular to velocity, but confined to 2D
        lift_direction = [-particles(i).velocity(2); particles(i).velocity(1)];
        lift_direction = lift_direction / (norm(lift_direction) + eps); % Normalize
        
        % Assign lift force
        particles(i).lift_force = lift_magnitude * lift_direction;

        % Update velocity using lift force
        particles(i).velocity = particles(i).velocity + timeStep * (particles(i).lift_force / m);

        % Update position
        particles(i).position = particles(i).position + particles(i).velocity * timeStep;

        % Check for boundary constraints
        for dim = 1:2 % Check both x and y dimensions
            if particles(i).position(dim) < bounds(dim, 1) || particles(i).position(dim) > bounds(dim, 2)
                particles(i).position(dim) = max(min(particles(i).position(dim), bounds(dim, 2)), bounds(dim, 1));
                particles(i).velocity(dim) = -cr * particles(i).velocity(dim); % Inelastic bounce
            end
        end

        % Draw pentagon
        verticesWorld = particles(i).vertices + particles(i).position;
        fill(verticesWorld(1, :), verticesWorld(2, :), 'b', 'FaceAlpha', 0.5);
        
        % Map particle area to grid
        x_min = particles(i).position(1) - pent_edge;
        x_max = particles(i).position(1) + pent_edge;
        y_min = particles(i).position(2) - pent_edge;
        y_max = particles(i).position(2) + pent_edge;
        
        % Find grid cells that the particle overlaps with
        x_cells = floor((x_min - bounds(1, 1)) / dx) + 1 : ceil((x_max - bounds(1, 1)) / dx);
        y_cells = floor((y_min - bounds(2, 1)) / dy) + 1 : ceil((y_max - bounds(2, 1)) / dy);
        
        % Ensure indices are within bounds
        x_cells = max(1, min(Nx, x_cells));
        y_cells = max(1, min(Ny, y_cells));
        
        % Add particle area to the grid cells
        for x_idx = x_cells
            for y_idx = y_cells
                area_fraction_grid(x_idx, y_idx) = area_fraction_grid(x_idx, y_idx) + A_pentagon;
            end
        end
    end
    
    % Normalize area fraction grid by cell area
    area_fraction_grid = area_fraction_grid / (dx * dy);
    
    % Visualize area fraction grid (optional)
    imagesc(bounds(1, :) + [0, Lx], bounds(2, :) + [0, Ly], area_fraction_grid');
    colorbar;
    title('Area Fraction Grid');
    
    % Handle Collisions and Update Velocities
    for i = 1:N_particles
        for j = i+1:N_particles
            if norm(particles(i).position - particles(j).position) < 2.0 * pent_edge
                % Compute collision normal
                normal = (particles(i).position - particles(j).position);
                normal = normal / norm(normal);
                
                % Compute relative velocity along the normal
                v_rel = dot(particles(i).velocity - particles(j).velocity, normal);
                
                if v_rel < 0 % Only resolve if approaching
                    % Compute impulse due to collision (in 2D)
                    J_c = -(1 + cr) * v_rel / (1/m + 1/m);
                    
                    % Apply impulse only in 2D
                    particles(i).velocity = particles(i).velocity + (J_c / m) * normal;
                    particles(j).velocity = particles(j).velocity - (J_c / m) * normal;
                end
            end
        end
    end
    
    drawnow;
    pause(0.01);
end