clear;
close;
clc;
format short g

%% Problem setup

n = 3;  % Number of customer vertices
s = 0;  % Number of ACS vertices
vehicles = 3;   % Number of vehicles
nf = 10; 
kilometers_per_second = 15;
T_max = 1000;
Q = 1161;
r = 1;

% syms v0
% syms v [1 n+s+(n*s)];
% v = [(1:n+s+nf*s)', randi(1000, n+s+nf*s, 2)];
rng(2);
v = [(1:n+s+nf*s)', rand(n+s+nf*s,1).*1000, rand(n+s+nf*s,1).*1000];
% v = [(1:n+s+nf*s)', [1:n+s+nf*s]', mod(linspace(0,100,n+s+nf*s)', 7)];
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
% F_0(size(F_0,1)+1, :) = [size(F_0,1), 0, 0]; % Add copy of depot
V = union(I,F, 'rows');
V = union(v0,V, 'rows');
V_prime = union(V, v(n+s+1:end,:), 'rows');
for i = 1:nf
    V_prime(size(V_prime,1)+1, :) = [size(V_prime,1), 0, 0]; % Add copy of depot
    F_0(end+1, :) = V_prime(end, :);
end
E = [0 0];
for i = 1:size(V_prime(:,1), 1)
    for j = 1:size(V_prime(:,1), 1)
        if i == 1 && j == 1
            continue
        end
        E(end+1, :) = [V_prime(i, 1), V_prime(j, 1)];
    end
end
% E = nchoosek(V_prime(:,1), 2); % Edges between each vertex

% Calculate the distance between each vertex
edges = E;
for i = 1:size(edges, 1)
    edge = edges(i,:);
    vertex_i = V_prime(edge(1)+1,2:3);
    vertex_j = V_prime(edge(2)+1,2:3);
    d(i,:) = pdist([vertex_i; vertex_j], 'euclidean');
end

% Time to cover distance
t = d ./ kilometers_per_second;

% Service time. Constant for everything
p = ones(size(edges,1), 1);

E = [E d t p];  % Associate distance, cost, and travel time between each edge


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
    "y_j", size(V_prime,1) + 1,...
    'Type',"continuous",...
    'LowerBound', 0,...
    "UpperBound",Q);


% Time variable specifying the time of arrival of a vehicle at
% vertex j, initialized to zero upon departure from the depot
tau = optimvar(...
    "tau", size(V_prime,1)+1,...
    'Type',"continuous",...
    'LowerBound', 0,...
    "UpperBound",T_max);

% Formulate problem
minprob = optimproblem;
% Objective Function (1)
cost = sum(x .* E(:,3));
minprob.Objective = cost;

% Constraint (2)
eqn2 = optimconstr(size(I,1));
idx = 1;
for i = I(:,1)'
    j = setdiff(V_prime(:,1), i);
    % x((ismember(E(:,1),i) & ismember(E(:,2),j))).show
    eqn2(idx) = sum(x(ismember(E(:,1),i) & ismember(E(:,2),j))) == 1;
    idx = idx + 1;
end
minprob.Constraints.constraints2 = eqn2;
% show(eqn2)

% Constraint (3)
eqn3 = optimconstr(size(F_0,1));
idx = 1;
for i = F_0(:,1)'
    j = setdiff(V_prime(:,1), i);
    eqn3(idx) = sum(x(ismember(E(:,1),i) & ismember(E(:,2),j))) <= 1;
    idx = idx + 1;
end
minprob.Constraints.constraints3 = eqn3;
% show(eqn3)

% Constraint (4)
eqn4 = optimconstr(size(V_prime,1));
idx = 1;
for j = V_prime(:,1)'
    i = V_prime(:,1);
    i = setdiff(i, j);
    eqn4(idx) = sum(x(ismember(E(:,1),j) & ismember(E(:,2),i))) ...
        - sum(x(ismember(E(:,1),i) & ismember(E(:,2),j))) == 0;
    idx = idx + 1;
end
minprob.Constraints.constraints4 = eqn4;
% show(eqn4)

% Constraint eq (5)
m = vehicles;
j = setdiff(V_prime(:,1), 0);
eqn5 = sum(x(ismember(E(:,1),0) & ismember(E(:,2),j))) <= m;
minprob.Constraints.constraints5 = eqn5;
% show(eqn5)

% Constraint eq(6)
j = setdiff(V_prime(:,1), 0);
eqn6 = sum(x(ismember(E(:,1),j) & ismember(E(:,2),0))) <= m;
minprob.Constraints.constraints6 = eqn6;
% show(eqn6)

% Constraint eq(7)
% tau variable indices are 1 more than the vertex they are associated with
eqn7 = optimconstr(size(V_prime,1) * 2 - 1);
idx = 1;
for i = V_prime(:,1)'
    for j = setdiff(V_prime(:,1), [0 i])'
        E_row = E(ismember(E(:,1),i) & ismember(E(:,2),j), :);
        eqn7(idx) = tau(j+1) >= tau(i+1) + ...
            (E_row(4) - E_row(5)) - ...
            T_max*(1 - x(ismember(E(:,1),i) & ismember(E(:,2),j)));
        idx = idx + 1;
    end
end
minprob.Constraints.constraints7 = eqn7;
% show(eqn7)

% Constraint eq(8)
% The first tau decision variable is for the depot
eqn8 = tau(1) <= T_max;
minprob.Constraints.constraints8 = eqn8;
% show(eqn8)

% Constraint eq(9)
eqn9_1 = optimconstr(size(V_prime, 1)-1);
eqn9_2 = optimconstr(size(V_prime, 1)-1);
idx = 1;
for j = setdiff(V_prime(:,1), 0)'
    eqn9_1(idx) = E(ismember(E(:,1),0) & ismember(E(:,2),j), 4) <= tau(j+1);
    eqn9_2(idx) = tau(j+1) <= T_max - ...
        (E(ismember(E(:,1),j) & ismember(E(:,2),0), 4) + ...
        E(ismember(E(:,1),j) & ismember(E(:,2),j), 5));
    idx = idx + 1;
end
minprob.Constraints.constraints9_1 = eqn9_1;
minprob.Constraints.constraints9_2 = eqn9_2;
% show(eqn9_1);
% show(eqn9_2);

% Constraint eq(10)
% y variable indices are 1 more than the vertex they are associated with
eqn10 = optimconstr(size(V_prime, 1) * size(I, 1));
idx = 1;
for i = V_prime(:,1)'
    for j = setdiff(I(:,1), i)'
        E_row = E(ismember(E(:,1),i) & ismember(E(:,2),j), :);
        eqn10(idx) = y(j+1) <= y(i+1) - ...
            r * E_row(3) * x(ismember(E(:,1),i) & ismember(E(:,2),j)) + ...
            Q * (1 - x(ismember(E(:,1),i) & ismember(E(:,2),j)));
        idx = idx + 1;
    end
end
minprob.Constraints.constraint10 = eqn10;
% show(eqn10);

% Constraint eq(11)
eqn11 = optimconstr(size(F_0, 1));
idx = 1;
for j = F_0(:,1)'
    eqn11(idx) = y(j+1) == Q;
    idx = idx + 1;
end
minprob.Constraints.constraint11 = eqn11;
% show(eqn11)

% Constraint eq(12)
idx = 1;
eqn12 = optimconstr(size(I,1)*size(F_prime,1));
for j = I(:,1)'
    for l = F_prime(:,1)'
        E_row_jl = E(ismember(E(:,1),j) & ismember(E(:,2),l), :);
        E_row_j0 = E(ismember(E(:,1),j) & ismember(E(:,2),0), :);
        E_row_l0 = E(ismember(E(:,1),l) & ismember(E(:,2),0), :);
        eqn12(idx) = y(j+1) >= min(r * E_row_j0(3),...
            r * (E_row_jl(3) + E_row_l0(3)));
        idx = idx + 1;
    end
end
minprob.Constraints.constraint12 = eqn12;
% show(eqn12);

T_max_delta = T_max;
while true
    % Solve
    [sol, fval] = solve(minprob);
    %[E(:,1:2), sol.x]
    sol.x = round(sol.x);
    tours = E(find(sol.x),1:2); % Filters to find edges travelled
    T_max_delta = T_max_delta/2;
    if T_max_delta <= 1 && ~isempty(tours)
        break
    elseif isempty(tours)
        T_max = T_max + T_max_delta;
    else
        T_max = T_max - T_max_delta;
    end
    tau.UpperBound = T_max;
 end
    
if isfield(sol,'tau')
    sol.tau
end
V_prime
T_max

figure;
hold on;
scatter(V_prime(1,2),V_prime(1,3),100,'k','filled');
scatter(V_prime(2:n+1,2),V_prime(2:n+1,3),100,'b','filled');
if s ~= 0
    scatter(V_prime(n+2:n+s+nf,2),V_prime(n+2:n+s+nf,3),100,'o','filled');
end
for i = tours'
    temp = [V_prime(ismember(V_prime(:,1),i(1)),2:3)', V_prime(ismember(V_prime(:,1),i(2)),2:3)'];
    plot(temp(1,:),temp(2,:),'k')
end
if s ~= 0
    legend("Depot","Customers","ACS");
else
    legend("Depot", "Customers");
end
axis([-10, 1000, -10, 1000]);