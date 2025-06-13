% Programmer Name: Justin Chen
% The following program will train a neural network with a set sample point
% and display the error from the actual data.

clear all; close all; clc;

% True trains the NN based on criterion below
% False load a NN trained on one of the two loaded examples
fromScratch = true;

% Select Sample Functions
% example = 'PaperLinear'; % f = @(X1, X2) cos(X1)+sin(X2);
example = '6_2linear'; % f = @(X1, X2) 3*X1;

% Define function you would like to train
% if (example == 'Papersincos')
    % f = @(X1, X2) cos(X1)+sin(X2)+2*X2;
    % h_layers = [20,10,10];      % Defining layers for NN
% else
    f = @(X1, X2) 3*X1+2*X2;
    h_layers = [4 5 4];         % Defining layers for NN
% end

if(fromScratch)
    % Training Size
    N = 100; % Larger N typically means longer training time
    ep = 1000;    % Number of Epochs
    
    % Training Paramters
    x1_train = linspace(-pi,pi,N);
    x2_train = linspace(-pi,pi,N);
    % x1_train = [linspace(-pi,0,N/2) linspace(0,pi,N/20)];
    % x2_train = [linspace(-pi,0,N/20) linspace(0,pi,N/2)];
    
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
    Y_train = f(X1_train,X2_train); % This data is the actual function mapped
    Y_validate = f(X1_test,X2_test);
    
    
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
               regressionLayer];
         
    options = trainingOptions("sgdm",'Momentum', momentum, MaxEpochs = ep,Plots='training-progress');
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
    load(append(name, ".mat"))
end

%% Larger training data that was not able to be plotted probably due to hardware limitations
% Trained on N = 200 and epoch of 200, 16 mins of training time

% load(200x200sincos_20_10_10.mat)

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

%% Plot neural network output space 

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

fprintf('Zonotope model: ')
toc

tempZono = NN.Z(NN.dimKeys);

subplot(1,4,4)
plot(tempZono,'r',0.8);
grid on;
title('Hybrid Zonotope')
toc
% %% Plot Hybrid Zonotope
% 
% % temp = zono([1 0; 0 1; 1 2],[0;0;0]);
% M = eye(3); M(3,2) = -2;
% tempZono = NN.Z;
% 
% tempM = M*NN.Z;
% subplot(1,4,4)
% plot(tempM,'r',0.8);
% grid on;
% title('Hybrid Zonotope')
% toc
%%
save(append(name, ".mat"))