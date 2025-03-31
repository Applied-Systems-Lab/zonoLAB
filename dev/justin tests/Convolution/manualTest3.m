clear
close all
clc

%%
load("GOODDONOTTOUCH_reLUMnistwConvPool(0.99)1111111.mat")

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
j=1;
% for j = 1:2163
    testPoint = reshape(XTest(:,:,1,j),[28 28])-net1.Layers(1).Mean;
    %% For the 7
    weights = net1.Layers(2,1).Weights;
    bias = net1.Layers(2,1).Bias;
    
    % figure
    % hold on
    
    % returns a 28 by 28 matrix
    for i = 1: length(net1.Layers(2).Weights)
        % subplot(10,10,i)
        % forward{i} = conv2(m7,reshape(weights(:,:,1,i),[5,5]),'valid') + bias(i);
        forward{i} = conv2(testPoint,reshape(rot90(weights(:,:,1,i),2),[5,5]),'same')+bias(i);
        % forward{i} = conv2(m7,reshape(rot90(weights(:,:,1,i),2),[5,5]),'valid') + bias(i);
        % imshow(forward{i})
        % title(i)
    end
    

    % Reshaping Input
    x=[];
    for i = 1:length(net1.Layers(2).Weights)
        x = [x ;reshape(forward{i},[1 size(forward{i},1)*size(forward{i},2)])];
    end

    %% Average Pooling Layer
    poolLength = net1.Layers(3).PoolSize(1);
    poolHeight = net1.Layers(3).PoolSize(2);
    strideVert = net1.Layers(3).Stride(1);
    strideHorz = net1.Layers(3).Stride(2);


    % Number of cells averaged
    pool = 1/(poolLength*poolHeight);

    k = 1;
    poolMatrix = pool*ones(net1.Layers(3).PoolSize);

    for j = 1: strideVert : size(forward{1},1)
        for i = 1: strideHorz :size(forward{1},2)
            zeroMatrix = zeros(size(forward{i}));
            zeroMatrix(j:j+poolHeight-1,i:i+poolLength-1) = zeroMatrix(j:j+poolHeight-1,i:i+poolLength-1)+poolMatrix;

            poolingShape(k,:) = reshape(zeroMatrix',[size(forward{i},1)*size(forward{i},2) 1])';
            k = k + 1;
        end
    end

    %%
    x = poolingShape*x';    % calcs fine
    %%
    % First Relu Layer (57600*50)
    % W * x + B
    
    %% Relu layer 4-5
    weights = net1.Layers(4,1).Weights;
    x = reshape(x,[size(weights,2) 1]);
    bias = net1.Layers(4,1).Bias;
    
    % linear mapping
    forwardRelu1 = weights*x+bias;
    
    % relu activation function
    activeRelu1 = max(0,forwardRelu1);
    
    %% Relu layer 6-7
    weights = net1.Layers(6,1).Weights;
    bias = net1.Layers(6,1).Bias;
    
    % linear mapping
    forwardRelu2 = weights*activeRelu1 + bias;
    
    % relu activation function
    activeRelu2 = max(0,forwardRelu2);
    
    %% Relu layer 8-9
    weights = net1.Layers(8,1).Weights;
    bias = net1.Layers(8,1).Bias;
    
    % linear mapping
    forwardRelu3 = weights*activeRelu2 + bias;
    
    % relu activation function
    activeRelu3 = max(0,forwardRelu3);
    
    %% Relu layer 10
    weights = net1.Layers(10,1).Weights;
    bias = net1.Layers(10,1).Bias;
    
    % linear mapping
    lastLayer(j) = weights*activeRelu3 + bias;
% end

%% Similarity Test