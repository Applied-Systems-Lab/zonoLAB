% clear
% clc
% close all


%% Getting online data
unzip("DigitsData.zip")

imds = imageDatastore("DigitsData", ...
    IncludeSubfolders=true, ...
    LabelSource="foldernames");

classNames = categories(imds.Labels);

%% Partition into different data sets

[imdsTrain,imdsValidation,imdsTest] = splitEachLabel(imds,0.7,0.15,0.15,"randomized");

%% training options
% type out

options = trainingOptions("sgdm", ...
    MaxEpochs=6, ...
    ValidationData=imdsValidation, ...
    ValidationFrequency=30, ...
    Plots="training-progress", ...
    Metrics="accuracy", ...
    Verbose=false);

%% Building the netowrk manually
inputSize = [28 28 1];
numClasses = 10;

filterSize = 3;
numFilters = 100;

layers = [
    imageInputLayer(inputSize)

    convolution2dLayer(filterSize, numFilters)
    batchNormalizationLayer
    reluLayer

    fullyConnectedLayer(numClasses)
    softmaxLayer
];

net1 = trainnet(imdsTrain,layers,"mse",options);

% If we define using a designer
% net1 = trainnet(imdsTrain,layers,"crossentropy",options);

%% Network testing
accuracy = testnet(net1,imdsValidation,"accuracy");

%% Predictions

scores = minibatchpredict(net1,imdsValidation);
YValidation = scores2label(scores,classNames);

%% Visualization

numValidationObservations = numel(imdsValidation.Files);
idx = randi(numValidationObservations,9,1);

figure
tiledlayout("flow")
for i = 1:9
    nexttile
    img = readimage(imdsValidation,idx(i));
    imshow(img)
    title("Predicted Class: " + string(YValidation(idx(i))))
end

%% Convo NN display
figure
tiledlayout("flow")

for i = 1:32
    % temp = net1.Layers(2,1).Weights(:,:,1,i);
    nexttile
    imagesc(net1.Layers(2,1).Weights(:,:,1,i))
    title("Node: " + i)
end
colormap(gray(256))

% %% Convo NN display Shifted by lowest
% figure
% tiledlayout("flow")
% title('Shifted by lowest')
% for i = 1:32
%     temp = net1.Layers(2,1).Weights(:,:,1,i);
%     nexttile
%     if(~isempty(abs(min(temp(temp<0)))))
%         imshow(net1.Layers(2,1).Weights(:,:,1,i)+abs(min(temp(temp<0))))
%     end
%     imshow(temp)
%     title("Node: " + i)
% end