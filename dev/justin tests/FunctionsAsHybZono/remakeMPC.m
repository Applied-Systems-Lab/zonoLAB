% Programmer Name: Justin Chen
% The following code is a updated version of the 2022-11-29 version of the
% MPC build hybzono command file form the L4DC papers.

% List of Changes
% Created a ReluNN in @memZono
% Created a new function 'relabel' that relabels factors, dims, and cons
% relabel function is mainly used so re can repeatedly use the same NN
% fucntion
% 
% Current Issue
% There is an overplot when you only plot out the X_all which will display
% the individual facets within each step in the X{i} evolutions. 

clear; close all;
clc;


%% Load training data
% Imports Weights, Ws and Biases, bs and a net for the NN
load MPC_NN_2_8_4_1.mat

%% System definitions

% Double integrator sytstem discretized with a sampling time of 1 second
A = [ 1 1;
      0 1];
B = [ 0.5;
      1];

n = size(A,1);  % size of A *probably assuming square

N = 6; % number of steps

%% Evolving the system
X = {};
U = {};             % Inputs
U_ = {};            % For intermediate calculation
% NN = {};

% Redefined initial set from original
% hybZono parameters are Zh = hybZono(Gc,Gb,c,Ac,Ab,b);
% Original Problem Statement has c as first parameter and Gc as second
% parameter

X{1} = hybZono(0.25*eye(2),[],[2.75;0],[],[],[]); % initial set
X{1} = memZono(X{1},'X_1');
X_all = X{1};

%% Creating the entire MPC policy as a memZono
centerMPC = [1 ; 0.5];
genMPC = [2, 0 ; 0, 1.5];

Xmpc= hybZono(genMPC, [], centerMPC, [], [], []);
Xmpc = memZono(Xmpc,'Xmpc');

%% Application of Relu NN function
% Usage of the reluNN function within the @memZono folder
[NNmpc] = reluNN(Xmpc, Ws, bs, 1000);
NN = NNmpc;

%% Time step Evolutions
for i = 1: N-1     % iterates through the defined steps
    i   % tracks the current iteration
    
    % New implementation of the relabel function that will work around the
    % reusing of a function
    % Found within the @memZono folder
    NN = NN.relabel(sprintf('_%d',i));

    % Reordering the dimision keys so they line up with U?
    NN.dimKeys = {X{i}.dimKeys{1}, X{i}.dimKeys{2}, sprintf('U_%d',i)};

    % Intersections current input X{i} with NN to get next U value
    U{i} = NN.merge(X{i}, sprintf('terminal_cons_%d',i));

    % After successfully merging we no longer need the information provided
    % by X{i}'s dim keys
    % "Simplification"
    U_{i} = U{i}({sprintf('U_%d',i)});

    % Calculates the next X{i}
    newDims = {sprintf('X_%d_1', i + 1),sprintf('X_%d_2',i + 1)};
    X{i + 1} = X{i}.transform(U_{i}.transform([],B,{},newDims),A,{},newDims);  % Transform has affine A*x + B

    % Stacking  
    X_all = X_all.merge(U_{i});
    X_all = X_all.merge(X{i+1});
end

%% Plotting All Figures
figure("name","All Plots")

plotTitle = ["X{i}","X(all)","All Plot"];

for i = 1: 3
    title(subplot(1,3,i),plotTitle(i));
    xlabel('x1');
    ylabel('x2');
    zlabel('u');
end

% Plots the basic X states
% Plots the facets that are present within the X states
% Plots all of it over the MPC policy to demonstrate the relation.
for i = 1: N-1     % iterates through the defined steps
    subplot(1,3,1)
    plot(X{i},{sprintf('X_%d_1', i),sprintf('X_%d_2',i)}, 'g', 1);

    subplot(1,3,2)
    plot(X_all,{sprintf('X_%d_1', i),sprintf('X_%d_2',i)}, 'b', 0.2);   % <- Plot gives the wrong color for the first facet?

    subplot(1,3,3)
    plot(U{i},{sprintf('X_%d_1', i),sprintf('X_%d_2',i), sprintf('U_%d',i)}, 'g', 1);
    plot(X{i},{sprintf('X_%d_1', i),sprintf('X_%d_2',i)}, 'g', 1);
    plot(X_all,{sprintf('X_%d_1', i),sprintf('X_%d_2',i)}, 'b', 0.5);
end

% Usage of .Z will need to be modified if the ordering of the dim keys are
% different.
subplot(1,3,3)
plot(NNmpc.Z, 'r', 0.8);