clc; clear; close all;

% Simulation Parameters
numParticles = 10;       % Number of pentagons
domainSize = 10;         % Simulation domain size
dt = 0.01;              % Time step
totalTime = 5;          % Total simulation time
fluid_on = true;        % Toggle fluid effect

% Physical Properties
mass = 1;               % Mass of each particle
sideLength = 1;         % Side length of pentagons
radius = sideLength / (2 * sin(pi/5)); % Circumradius of pentagon
C_L = 0.5;             % Lift coefficient (only for fluid case)
rho_fluid = 1.0;       % Fluid density
U_fluid = 2;           % Fluid velocity in y-direction

% Initialize Particles
positions = domainSize * rand(numParticles, 2); % Random positions
velocities = 2 * (rand(numParticles, 2) - 0.5); % Random initial velocities

% Ensure No Initial Overlap
for i = 1:numParticles
    for j = 1:i-1
        while norm(positions(i, :) - positions(j, :)) < 2 * radius
            positions(i, :) = domainSize * rand(1, 2);
        end
    end
end

% Simulation Loop
figure;
hold on;
axis equal;
xlim([0 domainSize]);
ylim([0 domainSize]);

for t = 0:dt:totalTime
    clf;
    hold on;
    xlim([0 domainSize]);
    ylim([0 domainSize]);

    % Compute Lift Force if Fluid is Present
    if fluid_on
        liftForces = C_L * rho_fluid * U_fluid^2 * ones(numParticles, 2);
        liftForces(:, 2) = liftForces(:, 2); % Lift in y-direction
    else
        liftForces = zeros(numParticles, 2); % No fluid force
    end

    % Update Velocities and Positions
    for i = 1:numParticles
        % Apply Lift Force
        acceleration = liftForces(i, :) / mass;
        velocities(i, :) = velocities(i, :) + acceleration * dt;
        positions(i, :) = positions(i, :) + velocities(i, :) * dt;

        % Collision Handling (Elastic Collision between Pentagons)
        for j = 1:numParticles
            if i ~= j
                dist = norm(positions(i, :) - positions(j, :));
                if dist < 2 * radius
                    % Compute collision normal
                    normal = (positions(i, :) - positions(j, :)) / dist;
                    relativeVel = velocities(i, :) - velocities(j, :);
                    
                    % Calculate the impulse
                    v_rel = dot(relativeVel, normal);
                    if v_rel < 0  % Only resolve if particles are approaching
                        impulse = -(1 + cr) * v_rel / (2 * mass);  % Elastic collision formula
                        
                        % Apply the impulse to the velocities
                        velocities(i, :) = velocities(i, :) + impulse * normal;
                        velocities(j, :) = velocities(j, :) - impulse * normal;

                        % Reposition particles slightly to prevent sticking
                        overlap = 2 * radius - dist;
                        correction = 0.5 * overlap * normal;
                        positions(i, :) = positions(i, :) + correction;
                        positions(j, :) = positions(j, :) - correction;
                    end
                end
            end
        end

        % Boundary Conditions (Reflective Walls)
        for d = 1:2
            if positions(i, d) - radius < 0 || positions(i, d) + radius > domainSize
                velocities(i, d) = -velocities(i, d);
            end
        end

        % Plot Regular Pentagon (No Rotation)
        theta = linspace(0, 2*pi, 6); % Pentagon shape
        x = radius * cos(theta) + positions(i, 1);
        y = radius * sin(theta) + positions(i, 2);
        fill(x, y, 'r');
    end

    pause(0.01);
end
