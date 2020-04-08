clear;
clc;

%% Problem setup

n = 5;  % Number of customer vertices
s = 5;  % Number of ACS vertices

syms v_0
syms v [1 n+s];

% v_0 is depot vertex. 
% v_1 through v_n are customer vertices. 
% v_n+1 through v_s are ACS vertices
I = v(1:n);
F = v(n+1:n+s);
V = union(I,F);
V = union(v_0,V);
E = nchoosek(V, 2) % Edges between each vertex
syms c [1 size(E,1)]    % Cost of each edge
syms d [1 size(E,1)]    % Distance of each edge

t = rand(1, size(E,1))';
%% Problem Formulation
% Set up decision variables

% Decision variable to travel between vertices i and j
x = optimvar(...
    "x", size(E,1),...
    'Type', 'integer',...
    'LowerBound', 0,...
    "UpperBound", 1);

% Energy level variable specifying the remaining battery energy
% level upon arrival to vertex j. It is reset to Q at each
% charging station vertex i and the depot
y_j = optimvar(...
    "y_j", size(V,1),...
    'Type',"continuous",...
    'LowerBound', 0,...
    "UpperBound",10);


% Time variable specifying the time of arrival of a vehicle at
% vertex j, initialized to zero upon departure from the depot
tau_j = optimvar(...
    "tau_j", size(V,1),...
    'Type',"continuous",...
    'LowerBound', 0,...
    "UpperBound",1000);
% Formulate problem

minprob = optimproblem;
% Objective Function (1)
% 

cost = x .* t;
minprob.Objective = sum(cost);
% Constraint (2)
% 


% Solve

[sol, fval] = solve(minprob)