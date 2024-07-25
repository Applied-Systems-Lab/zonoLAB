% Programmer Name: Justin Chen
% The following program is will train a NN based on a given equation,
% theoretically increasingly higher dimension

clear; close all;
clc;

% True trains the NN based on criterion below
% False load a NN trained on one of the two loaded examples
fromScratch = true;

%% Create Training & Testing Data
example = 'hyperplane';
N = 100;                        % Sampling size
f = @(X1,X2,X3) 4*X1 + 2*X2 + 5*X3;   % equation
minVal= -3;
maxVal = 3;
numFeatures = 2;                % Number of inputs

if(fromScratch)
    % Training Paramters
    x1_train = linspace(minVal,maxVal,N);
    x2_train = linspace(minVal,maxVal,N);
    x3_train = linspace(minVal,maxVal,N);
    
    [X1_train, X2_train, X3_train] = meshgrid(x1_train,x2_train,x3_train);
    [n, m, o] = size(X1_train);
    input_train = [reshape(X1_train, n*m*o, 1), reshape(X2_train, n*m*o, 1), reshape(X3_train, n*m*o, 1)];
    
    % Testing Inputs
    x1_test = linspace(min(x1_train),max(x1_train),2*N);
    x2_test = x1_test;
    x3_test = x1_test;
    [X1_test, X2_test, X3_test] = meshgrid(x1_test,x2_test,x3_test);
    [p, q, r] = size(X1_test);
    input_test = [reshape(X1_test, p*q*r, 1), reshape(X2_test, p*q*r, 1), reshape(X3_test, p*q*r, 1)];
    
    % Training data and validation data
    Y_train = f(X1_train,X2_train, X3_train); % This data is the actual function mapped
    Y_validate = f(X1_test,X2_test,X3_test);
    
    
    %% Network Training
    output_train = reshape(Y_train,n*m*o,1); % reshape output data  
    
    % Within Training Options, between 0 and 1
    % Value determins how much influence the previous term has for the new
    % term. Defaults at 0.9.
    
    % At 0.95, NN will fail to build the linear function
    % Significance of the value?
    momentum = 0.85; 
    
    % setting up the actual network
    % variable definition to create varied layers for NN training
    layers = [featureInputLayer(3)];
    
    h_layers = [4,4,4];      % Defining layers for NN
    
    for i = 1: length(h_layers)
        layers = [ layers 
                    fullyConnectedLayer(h_layers(i))
                    reluLayer];
    end
    layers = [ layers 
               fullyConnectedLayer(1)
               regressionLayer];
    
         
    options = trainingOptions("sgdm",'Momentum', momentum, MaxEpochs = 100,Plots='training-progress');
    net = trainNetwork(input_train, output_train, layers, options);
    
    name = example;
    
    for i = 1: length(h_layers)
        name = append(name, sprintf('_%d',h_layers(i)))
    end
    
    save(append(name, ".mat"))
else
    name = example;

    for i = 1: length(h_layers)
        name = append(name, sprintf('_%d',h_layers(i)))
    end
    load(append(name, ".mat"),"net")
end

[x1_min, x1_max] = deal(double(min(x1_test)), double(max(x1_test)));
[x2_min, x2_max] = deal(double(min(x2_test)), double(max(x2_test)));
[x3_min, x3_max] = deal(double(min(x3_test)), double(max(x3_test)));
domain = [x1_min,x1_max,x2_min,x2_max,x3_min,x3_max];
Ws = [];
bs = [];

for i = 1: floor(length(layers)/2)
    Ws = [Ws {double(net.Layers(2*i).Weights)}];
    bs = [bs {double(net.Layers(2*i).Bias)}];
end

%% Construct Hybrid Zonotope
tic

% Construct input zonotope
cdomain = num2cell(domain);
[x1_min, x1_max, x2_min, x2_max,x3_min,x3_max] = cdomain{:};
g11 = (x1_max - x1_min)/2;
g22 = (x2_max - x2_min)/2;
g33 = (x3_max - x3_min)/2;
Gx = diag([g11, g22 g33]);
cx = zeros(3, 1);
X = hybZono(Gx, [], cx, [], [], []);
X = memZono(X,'X');
a = 1000;

[NN,Y] = reluNN(X,Ws,bs,a);