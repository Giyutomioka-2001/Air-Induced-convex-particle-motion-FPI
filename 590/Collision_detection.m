%% Hexagon Particle Interaction Simulation in MATLAB (No Overlap, No Collision, Stable Motion)

clc;
clear;
close all;

% Parameters
numParticles = 10;             % Number of particles
timeStep = 0.05;              % Time step for simulation
numSteps = 1000;              % Number of simulation steps
gravity = [0; -9.8];          % Gravitational acceleration [m/s^2]
bounds = [-10, 10; -10, 10];  % Simulation box bounds [x_min, x_max; y_min, y_max]
hexRadius = 1;                % Fixed size for all hexagons
mass = 1;                     % Uniform mass for all particles

% Generate non-overlapping initial positions
particles = struct();
positions = [];
for i = 1:numParticles
    while true
        pos = 8 * rand(2,1) - [4; 4]; % Random position
        if isempty(positions) || all(vecnorm(positions - pos, 2, 1) > 2.1*hexRadius)
            positions = [positions, pos];
            break;
        end
    end
    
    % Define hexagon vertices
    angle = linspace(0, 2*pi, 7);
    x = hexRadius * cos(angle);
    y = hexRadius * sin(angle);
    particles(i).vertices = [x; y];
    particles(i).position = pos;
    particles(i).velocity = [0; 0]; % Initial velocity is zero
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
                particles(i).velocity(dim) = 0; % Stop movement upon reaching boundary
            end
        end

        % Draw hexagon
        verticesWorld = particles(i).vertices + particles(i).position;
        fill(verticesWorld(1, :), verticesWorld(2, :), 'b', 'FaceAlpha', 0.5);
    end
    
    % Check for real-time interactions
    for i = 1:numParticles
        for j = i+1:numParticles
            verticesA = particles(i).vertices + particles(i).position;
            verticesB = particles(j).vertices + particles(j).position;
            
            % Edge-Edge Collision
            for k = 1:6
                for l = 1:6
                    if checkEdgeEdge(verticesA(:, k), verticesA(:, mod(k,6)+1), ...
                                     verticesB(:, l), verticesB(:, mod(l,6)+1))
                        disp(['Edge-Edge Collision Detected: Particle ', num2str(i), ' & ', num2str(j)]);
                    end
                end
            end
            
            % Edge-Vertex Collision
            for k = 1:6
                for l = 1:6
                    if checkEdgeVertex(verticesB(:, l), verticesA(:, k), verticesA(:, mod(k,6)+1), 0.1)
                        disp(['Edge-Vertex Collision Detected: Vertex of Particle ', num2str(j), ' on Edge of Particle ', num2str(i)]);
                    end
                end
            end
            
            % Vertex-Vertex Collision
            for k = 1:6
                for l = 1:6
                    if checkVertexVertex(verticesA(:, k), verticesB(:, l), 0.1)
                        disp(['Vertex-Vertex Collision Detected: Particle ', num2str(i), ' & ', num2str(j)]);
                    end
                end
            end
        end
    end
    
    drawnow;
    pause(0.01);
end

%% Function Definitions
function isIntersecting = checkEdgeEdge(p1, p2, q1, q2)
    function val = cross2D(a, b)
        val = a(1) * b(2) - a(2) * b(1);
    end
    d1 = cross2D(q2 - q1, p1 - q1);
    d2 = cross2D(q2 - q1, p2 - q1);
    d3 = cross2D(p2 - p1, q1 - p1);
    d4 = cross2D(p2 - p1, q2 - p1);
    isIntersecting = (d1 * d2 < 0) && (d3 * d4 < 0);
end

function isNear = checkEdgeVertex(vertex, edgeStart, edgeEnd, threshold)
    edgeVec = edgeEnd - edgeStart;
    vertexVec = vertex - edgeStart;
    projLength = dot(vertexVec, edgeVec) / norm(edgeVec);
    projPoint = edgeStart + (projLength / norm(edgeVec)) * edgeVec;
    distance = norm(projPoint - vertex);
    isNear = distance < threshold;
end

function isClose = checkVertexVertex(v1, v2, threshold)
    isClose = norm(v1 - v2) < threshold;
end
