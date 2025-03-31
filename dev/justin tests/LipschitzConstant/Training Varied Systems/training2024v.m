% Programmer Name: Justin Chen
% The following program is will train a NN based on a given equation,
% theoretically increasingly higher dimension

clear; close all;
clc;

% True trains the NN based on criterion below
% False load a NN trained on one of the two loaded examples
fromScratch = false;

%% Create Training & Testing Data
example = 'MIMO';

%% Equation Definition
% Currently working with 2 inputs is easiest in terms of indexing
% Adding rows to A increases outputs

% A = [1 3; 4 2];
A = [1 9; 8 2;5 5];

f = @(X1,X2) A*[X1;X2];   % equation

numInputs = size(A,2);
numOutputs = size(A,1);

h_layers = [8 4 3];     % Defining layers for NN

N = 100;                % Sampling size Larger means longer training

minVal= -5;
maxVal = 5;

%% Naming for Data to be generated
name = append(sprintf('%din-%dout-',numInputs,numOutputs), example);

for i = 1: length(h_layers)
    name = append(name, sprintf('_%d',h_layers(i)))
end

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
    Y_train = arrayfun(f,X1_train,X2_train,'UniformOutput',false); % This data is the actual function mapped
    Y_validate = arrayfun(f,X1_test,X2_test,'UniformOutput',false);
    
    
    %% Network Training
    output_train = horzcat(Y_train{:})'; % reshape output data  
    output_validate = horzcat(Y_validate{:})'; % reshape output data  
    
    % Within Training Options, between 0 and 1
    % Value determins how much influence the previous term has for the new
    % term. Defaults at 0.9.
    
    % At 0.95, NN will fail to build the linear function
    % Significance of the value?
    momentum = 0.95; 
    
    % setting up the actual network
    % variable definition to create varied layers for NN training
    layers = [featureInputLayer(numInputs)];
    
    for i = 1: length(h_layers)
        layers = [ layers 
                    fullyConnectedLayer(h_layers(i))
                    reluLayer];
    end
    layers = [ layers 
               fullyConnectedLayer(numOutputs)];
    
         
    options = trainingOptions("sgdm",'Momentum', momentum, MaxEpochs = 200,Plots='training-progress');
    net = trainnet(input_train, output_train, layers, "mae", options);
else
    load(append(name, ".mat"))
end

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

inDims = NN.dimKeys(1:numInputs);
outDims = NN.dimKeys(numInputs + 1 : numInputs + numOutputs);

save(append(name, ".mat"))

output_test = predict(net,input_test);

% Reshaping outputs

for i = 1: numOutputs
    Y_tested{i} = reshape(output_test(:,i),p,q);
    
    Y_trained{i} = reshape(output_train(:,i),n, m);
    
    Y_validated{i} = reshape(output_validate(:,i),p,q);

    outputNN{i} = NN.projection([inDims outDims(i)]);
end

%% Display f based on data points
% Similar to normal graphing methods
figure('Name',name)

for j = 1: numOutputs
    row = j - 1;
    % Acutal Function
    subplot(numOutputs,4,1+(row*4))
    surf(X1_train, X2_train, Y_trained{j}, 'EdgeColor', 'none');
    title(sprintf('Actual Function %d',j));

    % Neural Net Approximation
    subplot(numOutputs,4,2+(row*4))
    surf(X1_test, X2_test, Y_tested{j}, 'EdgeColor', 'none');
    title(sprintf('NN Approx Function %d',j));

    % Difference between the first and second plot
    subplot(numOutputs,4,3+(row*4))
    surf(X1_test,X2_test,Y_tested{j}-Y_validated{j}, 'EdgeColor', 'none');
    title(sprintf('Approximation Error %d',j));

    % Plot as Zonotope
    subplot(numOutputs,4,4+(row*4))
    plot(outputNN{j}.Z, 'r', '0.8');
    title(sprintf('Hybrid Zonotope Function %d',j));
    
    avgEr{j} = rmse(Y_tested{j},Y_validated{j},'all');
    fprintf('RMSE Function %d: %d\n', j, avgEr{j});
end