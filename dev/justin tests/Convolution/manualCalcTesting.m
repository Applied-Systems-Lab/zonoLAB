clear
close all
clc

%%
load("manualCalculation.mat")

%%
m7 = reshape(XTest(:,:,1,1),[28 28]);
m7l = XTestLabels(1);
m1 = reshape(XTest(:,:,1,2),[28 28]);
m1l = XTestLabels(2);

%% Step by step calculation through NN
newNet = dlnetwork;
layerProjection = [];

layers = [
    imageInputLayer([28 28 1],Mean = net1.Layers(1).Mean)]
newNet = addLayers(newNet,layers);
newNet = initialize(newNet)
newNet.predict(m7)

for i = 1: length(net1.Layers)-1
    layers = net1.Layers(i+1)
    newNet = addLayers(newNet,layers)
    newNet = connectLayers(newNet, newNet.Layers(i).Name, newNet.Layers(i+1).Name)
    newNet = initialize(newNet)
    layerProjection{i} = newNet.predict(m7)
    plot(newNet)
end

m7 = reshape(XTest(:,:,1,1),[28 28])-net1.Layers(1).Mean;
m1 = reshape(XTest(:,:,1,2),[28 28])-net1.Layers(1).Mean; %

%%
figure
imshow(m7)
figure
imshow(m1)
j=1
% for j = 1:2163
    testPoint = reshape(XTest(:,:,1,j),[28 28])-net1.Layers(1).Mean;
    %% For the 7
    weights = net1.Layers(2,1).Weights;
    bias = net1.Layers(2,1).Bias;
    
    % figure
    % hold on
    
    % returns a 24 by 24 matrix
    for i = 1: 100
        % subplot(10,10,i)
        % forward{i} = conv2(m7,reshape(weights(:,:,1,i),[5,5]),'valid') + bias(i);
        forward{i} = conv2(testPoint,reshape(rot90(weights(:,:,1,i),2),[5,5]),'valid')+bias(i);
        % forward{i} = conv2(m7,reshape(rot90(weights(:,:,1,i),2),[5,5]),'valid') + bias(i);
        % imshow(forward{i})
        % title(i)
    end
    
    %%
    % First Relu Layer (57600*50)
    % W * x + B
    x=[];
    for i = 1:100
        x = [x ;reshape(forward{i},[576,1])];
    end
    
    %% Relu layer 3-4
    weights = net1.Layers(3,1).Weights;
    bias = net1.Layers(3,1).Bias;
    
    % linear mapping
    forwardRelu1 = weights*x+bias;
    
    % relu activation function
    activeRelu1 = max(0,forwardRelu1);
    
    %% Relu layer 5-6
    weights = net1.Layers(5,1).Weights;
    bias = net1.Layers(5,1).Bias;
    
    % linear mapping
    forwardRelu2 = weights*activeRelu1 + bias;
    
    % relu activation function
    activeRelu2 = max(0,forwardRelu2);
    
    %% Relu layer 7-8
    weights = net1.Layers(7,1).Weights;
    bias = net1.Layers(7,1).Bias;
    
    % linear mapping
    forwardRelu3 = weights*activeRelu2 + bias;
    
    % relu activation function
    activeRelu3 = max(0,forwardRelu3);
    
    %% Relu layer 9
    weights = net1.Layers(9,1).Weights;
    bias = net1.Layers(9,1).Bias;
    
    % linear mapping
    lastLayer(j) = weights*activeRelu3 + bias;
    
    % relu activation function
    % classification(j,1) = lastLayer;
% end

%% Similarity Test

size(find(isapprox(lastLayer',net1.predict(XTest),'loose')==1))