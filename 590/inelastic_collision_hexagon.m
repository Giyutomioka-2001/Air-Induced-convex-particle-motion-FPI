clc;
clear;
close all;

% Parameters
numParticles = 10;             % Number of particles
timeStep = 0.05;               % Time step for simulation
numSteps = 1000;               % Number of simulation steps
gravity = [0; -9.8];           % Gravitational acceleration [m/s^2]
bounds = [-10, 10; -10, 10];   % Simulation box bounds [x_min, x_max; y_min, y_max]
pentagonRadius = 1;            % Fixed size for all pentagons
mass = 1;                      % Uniform mass for all particles
cr = 0.9;                      % Higher coefficient of restitution for more collisions

% Generate non-overlapping initial positions
particles = struct();
positions = [];
for i = 1:numParticles
    while true
        pos = 8 * rand(2,1) - [4; 4]; % Random position
        if isempty(positions) || all(vecnorm(positions - pos, 2, 1) > 1.0*pentagonRadius)
            positions = [positions, pos];
            break;
        end
    end
    
    % Define pentagon vertices (regular pentagon)
    angle = linspace(0, 2*pi, 6); % One extra point to close the pentagon
    x = pentagonRadius * cos(angle);
    y = pentagonRadius * sin(angle);
    particles(i).vertices = [x; y];
    particles(i).position = pos;
    particles(i).velocity = rand(2,1) * 2 - 1; % Random initial velocity to encourage multiple collisions
    particles(i).mass = mass;
end

% Simulation loop
figure;
for step = 1:numSteps
    clf;
    hold on;
    axis equal;
    xlim(bounds(1, :));
    ylim(bounds(2, :));
    
    for i = 1:numParticles
        % Update position under gravity
        particles(i).velocity = particles(i).velocity + gravity * timeStep;
        particles(i).position = particles(i).position + particles(i).velocity * timeStep;

        % Check for boundary constraints
        for dim = 1:2
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
    for i = 1:numParticles
        for j = i+1:numParticles
            if norm(particles(i).position - particles(j).position) < 2.1*pentagonRadius
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
