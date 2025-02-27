%% Pentagon Particle Interaction Simulation in MATLAB

clc;
clear;
close all;

% Parameters
numParticles = 5;           % Number of particles
timeStep = 0.01;              % Time step for simulation
numSteps = 1000;              % Number of simulation steps
gravity = [0; 0];          % Gravitational acceleration [m/s^2]
restitution = 0.8;            % Coefficient of restitution (elasticity)
bounds = [-10, 10; -10, 10]; % Simulation box bounds [x_min, x_max; y_min, y_max]

% Initialize particle properties
particles = struct();
for i = 1:numParticles
    % Generate solid pentagon particles
    numSides = 5; % Fixed to pentagon
    radius = rand() + 0.5; % Random size
    angle = linspace(0, 2*pi, numSides + 1); % Pentagon vertices
    x = radius * cos(angle(1:end-1));
    y = radius * sin(angle(1:end-1));
    particles(i).vertices = [x; y];
    particles(i).position = 5 * rand(2, 1) - [2.5; 2.5]; % Random position
    particles(i).velocity = rand(2, 1) - 0.5;            % Random velocity
    particles(i).mass = rand() + 0.5;                    % Random mass
    particles(i).volume = polyarea(x, y);                % Calculate volume of pentagon
end

% Simulation loop
figure;
for step = 1:numSteps
    % Clear current plot
    clf;
    hold on;
    axis equal;
    xlim(bounds(1, :));
    ylim(bounds(2, :));
    
    for i = 1:numParticles
        % Update position and velocity
        particles(i).velocity = particles(i).velocity + gravity * timeStep;
        particles(i).position = particles(i).position + particles(i).velocity * timeStep;

        % Check for collisions with walls
        for dim = 1:2
            if particles(i).position(dim) < bounds(dim, 1) || particles(i).position(dim) > bounds(dim, 2)
                particles(i).velocity(dim) = -restitution * particles(i).velocity(dim);
                particles(i).position(dim) = max(min(particles(i).position(dim), bounds(dim, 2)), bounds(dim, 1));
            end
        end

        % Draw particle
        verticesWorld = particles(i).vertices + particles(i).position;
        fill(verticesWorld(1, :), verticesWorld(2, :), 'b', 'FaceAlpha', 1.0); % Solid pentagon
    end

    % Check for particle-particle interactions (volume overlap)
    for i = 1:numParticles
        for j = i+1:numParticles
            % Check overlap using positions and volumes
            posDiff = particles(i).position - particles(j).position;
            dist = norm(posDiff);
            combinedRadius = max(particles(i).volume, particles(j).volume)^(1/2); % Approx radius from volume
            if dist < combinedRadius
                % Resolve overlap using elastic collision response
                v1 = particles(i).velocity;
                v2 = particles(j).velocity;
                m1 = particles(i).mass;
                m2 = particles(j).mass;

                particles(i).velocity = v1 - (2*m2 / (m1 + m2)) * (dot(v1 - v2, posDiff) / norm(posDiff)^2) * posDiff;
                particles(j).velocity = v2 - (2*m1 / (m1 + m2)) * (dot(v2 - v1, -posDiff) / norm(posDiff)^2) * -posDiff;

                % Slightly adjust positions to avoid sticking
                overlap = combinedRadius - dist;
                particles(i).position = particles(i).position + (overlap / 2) * (posDiff / dist);
                particles(j).position = particles(j).position - (overlap / 2) * (posDiff / dist);
            end
        end
    end

    % Update plot
    drawnow;
    pause(0.01);
end
