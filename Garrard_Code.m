clear;
clc;
format short g

%% Problem setup

n = 2;  % Number of customer vertices
s = 2;  % Number of ACS vertices
vehicles = 3;   % Number of vehicles
nf = 2;
meters_per_second = 15;
T_max = 125;

% syms v0
% syms v [1 n+s+(n*s)];
v = [(1:n+s+nf*s)', randi(1000, n+s+nf*s, 2)];
for i = n+1:size(v,1)-s
    v(s+i,2:end) = v(i,2:end);
end
v0 = [0 0 0];

% v_0 is depot vertex. 
% v_1 through v_n are customer vertices. 
% v_n+1 through v_s are ACS vertices
I = v(1:n,:);
I_0 = union(I, v0, 'rows');
F = v(n+1:n+s,:);
F_prime = union(F, v(n+s+1:end,:), 'rows');
F_0 = union(F_prime, v0, 'rows');
V = union(I,F, 'rows');
V = union(v0,V, 'rows');
V_prime = union(V, v(n+s+1:end,:), 'rows');
for i = size(V_prime,1)+1:size(V_prime,1)+nf
    V_prime(i,:) = [i-1, 0, 0]; % Add copies of depot
end
E = nchoosek(V_prime(:,1), 2); % Edges between each vertex

% Calculate the distance between each vertex
edges = nchoosek(V_prime(:,1), 2);
for i = 1:size(edges, 1)
    edge = edges(i,:);
    vertex_i = V_prime(edge(1)+1,2:3);
    vertex_j = V_prime(edge(2)+1,2:3);
    d(i,:) = pdist([vertex_i; vertex_j], 'euclidean');
end

% Time to cover distance
t = d ./ meters_per_second;

% Service time. Constant for everything
p = ones(size(edges,1), 1);

E = [E d t p];  % Associate distance, cost, and travel time between each edge

% Expand E so that each edge is bidirectional
edge_i = E(:,1);
edge_j = E(:,2);
E_temp = E;
E_temp(:, 1) = edge_j; 
E_temp(:, 2) = edge_i;
E = [E;E_temp];

%% Problem Formulation

% Set up decision variables

% Decision variable to travel between vertices i and j. Equivalent to eq
% (13)
x = optimvar(...
    "x", size(E,1),...
    'Type', 'integer',...
    'LowerBound', 0,...
    "UpperBound", 1);

% Energy level variable specifying the remaining battery energy
% level upon arrival to vertex j. It is reset to Q at each
% charging station vertex i and the depot
y = optimvar(...
    "y_j", size(V_prime,1),...
    'Type',"continuous",...
    'LowerBound', 0,...
    "UpperBound",10);


% Time variable specifying the time of arrival of a vehicle at
% vertex j, initialized to zero upon departure from the depot
tau = optimvar(...
    "tau", size(V_prime,1)+1,...
    'Type',"continuous",...
    'LowerBound', 0,...
    "UpperBound",1000);

% Formulate problem
minprob = optimproblem;
% Objective Function (1)
cost = x .* E(:,3);
minprob.Objective = sum(cost);

% Constraint (2)
eqn2 = optimconstr(size(I,1));
for i = 1:size(I,1)
    eqn2(i) = sum(x(ismember(E(:,1),I(i,1)) & ismember(E(:,2),V_prime(:,1)))) == 1;
end
minprob.Constraints.constraints2 = eqn2;
% show(eqn2)

% Constraint (3)
eqn3 = optimconstr(size(F_0,1));
for i = 1:size(F_0,1)
    eqn3(i) = sum(x(E(:,1) == F_0(i,1))) <= 1;
end
minprob.Constraints.constraints3 = eqn3;
% show(eqn3)

% Constraint (4)
eqn4 = optimconstr(size(V_prime,1));
for i = 1:size(V_prime,1)
    eqn4(i) = sum(x(E(:,1) == V_prime(i,1))) - sum(x(E(:,2) == V_prime(i,1))) == 0;
end
minprob.Constraints.constraints4 = eqn4;
% show(eqn4)

% Constraint eq (5)
m = vehicles;
eqn5 = sum(x(E(:,1) == 0)) <= m;
minprob.Constraints.constraints5 = eqn5;
% show(eqn5)

% Constraint eq(6)
eqn6 = sum(x(E(:,2) == 0)) <= m;
minprob.Constraints.constraints6 = eqn6;
% show(eqn6)

% Constraint eq(7).
V_prime_without_depot = setdiff(V_prime, v0, 'rows');
eqn7 = optimconstr(size(V_prime,1) * size(V_prime_without_depot,1));
for i = 1:size(V_prime,1)
    for j = 1:size(V_prime_without_depot,1)
        if j == i
            continue
        end
        % i-1
        % j-1
        E_row_idx = find(E(:,1) == i-1 & E(:,2) == j-1);
        E_row = E(E(:,1) == i-1 & E(:,2) == j-1, :);
        eqn7((i-1) * size(V_prime,1) + j) = tau(j) >= tau(i) ...
            + (E_row(4) - E_row(5)) * x(E_row_idx) ...
            - T_max * (1 - x(E_row_idx));
        % eqn7((i-1) * size(V_prime,1) + j).show;
    end
end
minprob.Constraints.constraints7 = eqn7;
% show(eqn7);

% Constraint eq(8)
eqn8 = tau(1) <= T_max;
minprob.Constraints.constraints8 = eqn8;
% show(eqn8)

% Constraint eq(9)
eqn9_1 = optimconstr(size(V_prime_without_depot, 1));
eqn9_2 = optimconstr(size(V_prime_without_depot, 1));
for j = 1:size(V_prime_without_depot, 1)
    E_row = E(E(:,1) == 0 & E(:,2) == j, :);
    eqn9_1(j) = E_row(4) <= tau(j+1);
    E_row = E(E(:,1) == j & E(:,2) == 0, :);
    eqn9_2(j) = tau(j+1) <= T_max - (E_row(4) + E_row(5));
    % eqn9_1(j).show
    % eqn9_2(j).show
end
minprob.Constraints.constraint9_1 = eqn9_1;
minprob.Constraints.constraint9_2 = eqn9_2;
% show(eqn9_1)
% show(eqn9_2)

% Solve
[sol, fval] = solve(minprob);
[E(:,1:2), sol.x]
sol.tau