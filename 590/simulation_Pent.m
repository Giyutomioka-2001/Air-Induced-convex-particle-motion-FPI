% MATLAB Code for Simulating Pentagon-Shaped Particles in a Circular Ring (2D)
% Particles in circular motion
clc; clear; close all;

%% Simulation Parameters
N_particles = 15; % Number of pentagon particles
radius = 20; % Radius of circular domain in cm
pent_edge = 1.0 ; % Edge length of the pentagon in cm
bounds = [-10, 10; -10, 10]; % Simulation box bounds [x_min, x_max; y_min, y_max]
mass = 1; % Uniform mass for all particles
cr = 0.9; % Higher coefficient of restitution for more collisions

Cl = 0.8; % Assumed moderate lift coefficient
V_air = 1.0; % Moderate air velocity
rho_air = 1.2; % Air density SI units

% Compute apothem
apothem = pent_edge / (2 * tan(pi/5)); % Apothem for a regular pentagon
A_pentagon = (5/2) * pent_edge * apothem; % Approx. area of pentagon

timeStep = 0.2; % Time step for simulation
numSteps = 2500; % Number of simulation steps

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
    particles(i).mass = mass;
    particles(i).lift_force = [0; 0]; % Initialize lift force
end

%% Create Video Writer
%videoFileName = 'pentagon_simulation.avi'; % Output video file name
%vidObj = VideoWriter(videoFileName, 'Uncompressed AVI'); % Create video object
%open(vidObj); % Open the video file

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
        if v_mag > 0
            velocity_unit = particles(i).velocity / v_mag;
            lift_direction = [-velocity_unit(2); velocity_unit(1)]; % Rotate velocity vector 90°
        else
            lift_direction = [0; 0]; % No lift if no motion
        end

        % Compute lift force magnitude
        lift_magnitude = Cl * (rho_air * V_air^2 / 2) * A_pentagon;
        particles(i).lift_force = [0; lift_magnitude]; % Lift only in +y direction

        % Update velocity using lift force (no rotational effects)
        particles(i).velocity = particles(i).velocity + timeStep * particles(i).lift_force / mass;

        % Update position due to velocity
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
    end
    
    % Check for multiple collisions and resolve them (only normal impulse)
    for i = 1:N_particles
        for j = i+1:N_particles
            if norm(particles(i).position - particles(j).position) < 2.1 * pent_edge
                % Compute collision normal
                normal = (particles(i).position - particles(j).position);
                normal = normal / norm(normal);
                
                % Relative velocity along the normal
                v_rel_normal = dot(particles(i).velocity - particles(j).velocity, normal);
                
                if v_rel_normal < 0 % Only resolve if approaching
                    % Compute impulse due to collision (normal impulse only)
                    J_c = -(1 + cr) * v_rel_normal / (1/mass + 1/mass);
                    
                    % Apply impulse without rotation
                    particles(i).velocity = particles(i).velocity + (J_c / mass) * normal;
                    particles(j).velocity = particles(j).velocity - (J_c / mass) * normal;
                end
            end
        end
    end
    
    drawnow;
    
    % Write current frame to video
    frame = getframe(gcf); % Capture current figure
    writeVideo(vidObj, frame); % Write frame to video file
    
    pause(0.01);
end

% Close the video writer
close(vidObj);

fprintf('Simulation recorded in %s\n', videoFileName);
