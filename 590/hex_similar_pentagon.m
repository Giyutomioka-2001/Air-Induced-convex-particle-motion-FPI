% MATLAB Code for Simulating Pentagon-Shaped Particles in a Circular Ring (2D)

clc; clear; close all;

%% Simulation Parameters
N_particles = 10; % Number of pentagon particles
radius = 10; % Radius of the circular domain
pentagon_size = 1; % Side length of the pentagon
bounds = [-10, 10; -10, 10]; % Simulation box bounds [x_min, x_max; y_min, y_max]
mass = 1; % Uniform mass for all particles
cr = 0.9; % Higher coefficient of restitution for more collisions

Cl = 0.8; % Assumed moderate lift coefficient
V_air = 0.5; % Moderate air velocity
rho_air = 1.2; % Air density SI units
A_pentagon = (5/2) * pentagon_size^2 ; % Approx. area of pentagon

timeStep = 0.5; % Time step for simulation
numSteps = 1000; % Number of simulation steps

%% Initialize Particles
particles = struct();
positions = [];
for i = 1:N_particles
    while true
        pos = 8 * rand(2,1) - [1;5]; % Random position
        if isempty(positions) || all(vecnorm(positions - pos, 2, 1) > 1.0*pentagon_size)
            positions = [positions, pos];
            break;
        end
    end
    
    % Define pentagon vertices (regular pentagon)
    angle = linspace(0, 2*pi, 6); % One extra point to close the pentagon
    x = pentagon_size * cos(angle);
    y = pentagon_size * sin(angle);
    particles(i).vertices = [x; y];
    particles(i).position = pos;
    particles(i).velocity = rand(2,1) * 2 - 1; % Random initial velocity
    particles(i).mass = mass;
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
    
    for i = 1:N_particles
        % Compute lift force strictly using the given formula
        v_mag = norm(particles(i).velocity);
        lift_magnitude = Cl * (rho_air * V_air^2 / 2) * A_pentagon;
        lift_direction = [-particles(i).velocity(2); particles(i).velocity(1)];
        lift_direction = lift_direction / norm(lift_direction + eps); % Perpendicular to velocity
        particles(i).lift_force = lift_magnitude * lift_direction;

        % Update velocity using lift force
        particles(i).velocity = particles(i).velocity + timeStep * particles(i).lift_force / mass;

        % Update position
        particles(i).position = particles(i).position + particles(i).velocity * timeStep;

        % Check for boundary constraints
        for dim = 1:1
            if particles(i).position(dim) < bounds(dim, 1) || particles(i).position(dim) > bounds(dim, 2)
                particles(i).position(dim) = max(min(particles(i).position(dim), bounds(dim, 2)), bounds(dim, 1));
                particles(i).velocity(dim) = -cr * particles(i).velocity(dim); % Inelastic bounce
            end
        end

        % Draw pentagon
        verticesWorld = particles(i).vertices + particles(i).position;
        fill(verticesWorld(1, :), verticesWorld(2, :), 'b', 'FaceAlpha', 0.5);
    end
    
    % Check for multiple collisions and resolve them
    for i = 1:N_particles
        for j = i+1:N_particles
            if norm(particles(i).position - particles(j).position) < 2.1*pentagon_size
                % Compute collision normal
                normal = (particles(i).position - particles(j).position);
                normal = normal / norm(normal);
                
                % Relative velocity along the normal
                v_rel = dot(particles(i).velocity - particles(j).velocity, normal);
                
                if v_rel < 0 % Only resolve if approaching
                    % Compute impulse magnitude for multiple collisions
                    impulse = -(1 + cr) * v_rel / (1/mass + 1/mass);
                    
                    % Update velocities with impulse
                    particles(i).velocity = particles(i).velocity + (impulse / mass) * normal;
                    particles(j).velocity = particles(j).velocity - (impulse / mass) * normal;
                end
            end
        end
    end
    
    drawnow;
    pause(0.01);
end
