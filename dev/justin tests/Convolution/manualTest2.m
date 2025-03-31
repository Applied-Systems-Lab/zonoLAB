clear
clc
close all

load('reLUMnist.mat')

%%
m7 = reshape(XTest(:,:,1,1),[28 28]);
m7l = XTestLabels(1);
m1 = reshape(XTest(:,:,1,2),[28 28]);
m1l = XTestLabels(2);

m7predict = net1.predict(XTest(:,:,1,1));
m1predict = net1.predict(XTest(:,:,1,2));

x = reshape(m7,[784 1])-net1.Layers(1).Mean;

for i = 1: length(2)
    x = reshape(XTest(:,:,1,i),[784 1])-net1.Layers(1).Mean;
    %% For the 7
    weights2 = net1.Layers(2,1).Weights;
    bias2 = net1.Layers(2,1).Bias;
    
    % linear mapping
    forwardRelu1 = weights2*x+bias2;
    
    % relu activation function
    activeRelu1 = max(0,forwardRelu1);
    
    %% Relu layer 5-6
    weights4 = net1.Layers(4,1).Weights;
    bias4 = net1.Layers(4,1).Bias;
    
    % linear mapping
    forwardRelu2 = weights4*activeRelu1 + bias4;
    
    % relu activation function
    activeRelu2 = max(0,forwardRelu2);
    
    %% Relu layer 7-8
    weights6 = net1.Layers(6,1).Weights;
    bias6 = net1.Layers(6,1).Bias;
    
    % linear mapping
    forwardRelu3 = weights6*activeRelu2 + bias6;
    
    % relu activation function
    activeRelu3 = max(0,forwardRelu3);
    
    %% Relu layer 9
    weights8 = net1.Layers(8,1).Weights;
    bias8 = net1.Layers(8,1).Bias;
    
    % linear mapping
    lastLayer(i,1) = weights8*activeRelu3 + bias8;
    
    % relu activation function
    % classification(j,1) = lastLayer;
    % end
end
%% Step by step calculation through NN
newNet = dlnetwork;
layerProjection = [];

newNet = addLayers(newNet,net1.Layers(1));
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