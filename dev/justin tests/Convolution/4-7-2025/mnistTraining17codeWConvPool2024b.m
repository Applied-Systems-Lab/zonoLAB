clear
clc
close all

%%
load('mnist.mat')

ones_sevens_idx = logical((training.labels==7) + (training.labels==1)+ (training.labels==4));
training.images = training.images(:,:,ones_sevens_idx);
training.labels = training.labels(ones_sevens_idx);
numObservations = size(training.images,3);

[idxTrain,idxValidation] = trainingPartitions(numObservations,[0.85 0.15]);

XValidation = reshape(training.images(:,:,idxValidation),28, 28, 1, []);
XTrain = reshape(training.images(:,:,idxTrain),28, 28, 1, []);
XTrainLabels = double(training.labels(idxTrain,1) == [1,4,7]);
XValidLabels = double(training.labels(idxValidation,1) == [1,4,7]);

ones_sevens_idx = logical((test.labels==7) + (test.labels==1)+ (test.labels==4));
test.images = test.images(:,:,ones_sevens_idx);
test.labels = test.labels(ones_sevens_idx);
numObservations = size(test.images,3);
XTest = reshape(test.images(:,:,:),28, 28, 1, []);
XTestLabels = double(test.labels(:,1) == [1,4,7]);

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
numClasses = 3;

filterSize = 5;
numFilters = 20;

layers = [
    imageInputLayer([28 28 1])

    % MUST ADD PADDING
    convolution2dLayer(filterSize, 10*numFilters,Padding=[filterSize-1 filterSize-1])      % should be full convolution

    %Need to add stride in order have no overlap
    averagePooling2dLayer(2,'Stride',2)

    fullyConnectedLayer(100)
    reluLayer

    fullyConnectedLayer(50)
    reluLayer

    fullyConnectedLayer(25)
    reluLayer

    fullyConnectedLayer(10)
    reluLayer

    fullyConnectedLayer(numClasses)
];

net1 = trainnet(XTrain, XTrainLabels,layers,"mae",options);

% If we define using a designer
% net1 = trainnet(imdsTrain,layers,"crossentropy",options);

%% Network testing
test_predictions = predict(net1,XTest);

%% Validation

% Max value per row, corresponding to highest confidence
classAboveConf = max(test_predictions,[],2) > 0.5; 

% Accuracy not accounting for treshold
classCorr = (XTestLabels - (test_predictions == max(test_predictions,[],2)));
accuracySoft = (1 - length(find(classCorr == -1))/size(classCorr,1))*100;

% Accuracy accounting for treshold
classBelowTreshold = XTestLabels - (~classAboveConf .* (test_predictions == max(test_predictions,[],2)));

% counts up all the rows that only have 0
% correct classification but not enough confidence
sum(all(classBelowTreshold == 0,2))

% finds all the -1s
% bad classification and not enough confidence
sum((classBelowTreshold == -1),'all')


classAboveTreshold = XTestLabels - (classAboveConf .* (test_predictions == max(test_predictions,[],2)));
% counts up all the rows that only have 0
% correct classification with enough confidence
sum(all(classAboveTreshold == 0,2))

% finds all the -1s
% bad classification with enough confidence
sum((classAboveTreshold == -1),'all')

%% Validation normalized
normTest_Predictions = normalize(test_predictions,"range");

% Max value per row, corresponding to highest confidence
classAboveConf = max(normTest_Predictions,[],2) > 0.5; 

% Accuracy not accounting for treshold
classCorr = (XTestLabels - (normTest_Predictions == max(normTest_Predictions,[],2)));
accuracySoft = (1 - length(find(classCorr == -1))/size(classCorr,1))*100;

% Accuracy accounting for treshold
classBelowTreshold = XTestLabels - (~classAboveConf .* (normTest_Predictions == max(normTest_Predictions,[],2)));

% counts up all the rows that only have 0
% correct classification but not enough confidence
sum(all(classBelowTreshold == 0,2))

% finds all the -1s
% bad classification and not enough confidence
sum((classBelowTreshold == -1),'all')


classAboveTreshold = XTestLabels - (classAboveConf .* (normTest_Predictions == max(normTest_Predictions,[],2)));
% counts up all the rows that only have 0
% correct classification with enough confidence
sum(all(classAboveTreshold == 0,2))

% finds all the -1s
% bad classification with enough confidence
sum((classAboveTreshold == -1),'all')