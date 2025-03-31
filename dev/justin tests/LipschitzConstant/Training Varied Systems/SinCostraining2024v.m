% Programmer Name: Justin Chen
% The following program is will train a NN based on a given equation,
% theoretically increasingly higher dimension

clear; close all;
clc;

% True trains the NN based on criterion below
% False load a NN trained on one of the two loaded examples
fromScratch = false;

%% Create Training & Testing Data
example = 'SinCosV2';
h_layers = [20,10,10];      % Defining layers for NN

N = 100;                        % Sampling size

%% Equation Definition
f = @(X1,X2) sin(X1)+cos(X2);   % equation
minVal= -5;
maxVal = 5;

if(fromScratch)
    % Training Paramters
    x1_train = linspace(minVal,maxVal,N);
    x2_train = linspace(minVal,maxVal,N);
    
    [X1_train, X2_train] = meshgrid(x1_train,x2_train);
    [n, m] = size(X1_train);
    input_train = [reshape(X1_train, n*m, 1), reshape(X2_train, n*m, 1)];
    
    % Testing Inputs
    x1_test = linspace(min(x1_train),max(x1_train),2*N);
    x2_test = x1_test;
    [X1_test, X2_test] = meshgrid(x1_test,x2_test);
    [p, q] = size(X1_test);
    input_test = [reshape(X1_test, p*q, 1), reshape(X2_test, p*q, 1)];
    
    % Training data and validation data
    Y_train = arrayfun(f,X1_train,X2_train,'UniformOutput',true); % This data is the actual function mapped
    Y_validate = arrayfun(f,X1_test,X2_test,'UniformOutput',true);  % For multiple output, uniform output false
    
    
    %% Network Training
    output_train = reshape(Y_train,n*m,1); % reshape output data  
    
    % Within Training Options, between 0 and 1
    % Value determins how much influence the previous term has for the new
    % term. Defaults at 0.9.
    
    % At 0.95, NN will fail to build the linear function
    % Significance of the value?
    momentum = 0.95; 
    
    % setting up the actual network
    % variable definition to create varied layers for NN training
    layers = [featureInputLayer(2)];
    
    for i = 1: length(h_layers)
        layers = [ layers 
                    fullyConnectedLayer(h_layers(i))
                    reluLayer];
    end
    layers = [ layers 
               fullyConnectedLayer(1)
               ];
    
         
    options = trainingOptions("sgdm",'Momentum', momentum, MaxEpochs = 200,Plots='training-progress');
    net = trainnet(input_train, output_train, layers, "mse", options);
    
    name = example;

    for i = 1: length(h_layers)
        name = append(name, sprintf('_%d',h_layers(i)))
    end
else
    name = example;

    for i = 1: length(h_layers)
        name = append(name, sprintf('_%d',h_layers(i)))
    end
    load(append(name, ".mat"))
end

%% Display f based on data points
% Similar to normal graphing methods
figure('Name',name)
subplot(1,4,1)
surf(X1_train, X2_train, Y_train, 'EdgeColor', 'none');
title('Actual Function')

%% Validate the Network
output_test = predict(net,input_test); 
Y_test = reshape(output_test,p,q);

% The approximation calculated from the training
subplot(1,4,2)
surf(X1_test,X2_test,Y_test, 'EdgeColor', 'none');
title('Neural Net Approx')

% Difference between the first and second plot
subplot(1,4,3)
surf(X1_test,X2_test,Y_test-Y_validate, 'EdgeColor', 'none');
title('Approximation Error')

%% Calculate Average Error
% Observe how bad the error is across the entire mapping
avgEr = rmse(Y_test,Y_validate,'all');
fprintf('RMSE: %d\n',avgEr);

%% Plot NN output space
[x1_min, x1_max] = deal(double(min(x1_test)), double(max(x1_test)));
[x2_min, x2_max] = deal(double(min(x2_test)), double(max(x2_test)));
domain = [x1_min,x1_max,x2_min,x2_max];
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
[x1_min, x1_max, x2_min, x2_max] = cdomain{:};
g11 = (x1_max - x1_min)/2;
g22 = (x2_max - x2_min)/2;
Gx = diag([g11, g22]);
cx = zeros(2, 1);
X = hybZono(Gx, [], cx, [], [], []);
X = memZono(X,'X');
a = 1000;

[NN,Y] = reluNN(X,Ws,bs,a);

%% Plot Hybrid Zonotope
subplot(1,4,4)
plot(NN.Z,'r',0.8);
grid on;
title('Hybrid Zonotope')
toc

save(append(name, ".mat"))