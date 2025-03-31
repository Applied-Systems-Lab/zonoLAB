clear
clc
close all


% %% Getting online data
% % unzip("DigitsData.zip")
% 
% imds = imageDatastore("DigitsData", ...
%     IncludeSubfolders=true, ...
%     LabelSource="foldernames");
% 
% classNames = categories(imds.Labels);
% 
% %% Partition into different data sets
% 
% [imdsTrain,imdsValidation,imdsTest] = splitEachLabel(imds,0.7,0.15,0.15,"randomized");

%%
load('mnist.mat')

ones_sevens_idx = logical((training.labels==7) + (training.labels==1));
training.images = training.images(:,:,ones_sevens_idx);
training.labels = training.labels(ones_sevens_idx);
numObservations = size(training.images,3);

[idxTrain,idxValidation] = trainingPartitions(numObservations,[0.85 0.15]);

XValidation = reshape(training.images(:,:,idxValidation),28, 28, 1, []);
XTrain = reshape(training.images(:,:,idxTrain),28, 28, 1, []);
XTrainLabels = training.labels(idxTrain,1);
XValidLabels = training.labels(idxValidation,1);

ones_sevens_idx = logical((test.labels==7) + (test.labels==1));
test.images = test.images(:,:,ones_sevens_idx);
test.labels = test.labels(ones_sevens_idx);
numObservations = size(test.images,3);
XTest = reshape(test.images(:,:,:),28, 28, 1, []);
XTestLabels = test.labels(:,1);

%% training options
% type out

options = trainingOptions("sgdm", ...
    MaxEpochs=10, ...
    ValidationData={XValidation,XValidLabels}, ...
    ValidationFrequency=30, ...
    Plots="training-progress", ...
    Metrics = "rmse", ...
    Verbose=false);

%% Building the netowrk manually
numClasses = 1;

filterSize = 3;
numFilters = 10;

layers = [
    imageInputLayer([28 28 1])

    % convolution2dLayer(filterSize, numFilters)

    fullyConnectedLayer(50)
    reluLayer

    fullyConnectedLayer(25)
    reluLayer

    fullyConnectedLayer(10)
    reluLayer

    fullyConnectedLayer(numClasses)
    % softmaxLayer
];

net1 = trainnet(XTrain, XTrainLabels,layers,"mse",options);

% If we define using a designer
% net1 = trainnet(imdsTrain,layers,"crossentropy",options);

%% Network testing
test_predictions = predict(net1,XTest);
prediction_error = abs( test_predictions - XTestLabels);              % Norm across each column

numClassError = length(find(prediction_error > 0.5));
accuracy = (1 - (numClassError/length(XTestLabels)))*100;

save('reLUMnist.mat')
%% Predictions

idx = randsample(length(XTestLabels),25);

figure
for i = 1 : size(idx,1)
    subplot(floor(sqrt(size(idx,1))),ceil(sqrt(size(idx,1))),i)      % Gets the cells of the worst 1
    imshow(test.images(:,:,idx(i)))
    if(prediction_error(idx(i)) > 0.5)
        title(sprintf('Wrong: %.3f', test_predictions(idx(i))))
    else
        title(sprintf('Right: %.3f', test_predictions(idx(i))))
    end
end

% Worst ones
worstidx = find(prediction_error > 0.5);

figure
for i = 1 : size(worstidx,1)
    subplot(floor(sqrt(size(worstidx,1))),ceil(sqrt(size(worstidx,1))),i)      % Gets the cells of the worst 1
    imshow(test.images(:,:,worstidx(i)))
    title(test_predictions(worstidx(i)))
end


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

%% Converting to hybzono



function [Z] = interval_hull_of_points(data,threshold)
    c = mean(data,2);
    M = max(data,[],2);
    m = min(data,[],2);
    s = std(data,0,2);
    if threshold == 'all'
        Gc = diag((M-m)/2);
    else
        Gc = diag(threshold*s);
    end
    Z = hybZono(Gc,[],c,[],[],[]);
    Z = memZono(Z,'Z');
end