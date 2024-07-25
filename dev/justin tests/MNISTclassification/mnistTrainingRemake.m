% Programmer Name: Justin Chen
% The following code is a updated version of the 2022-11-29 version of the
% mnist command file form the L4DC papers.

% List of Changes


clear; close all;
clc;

% True trains the mnist data
% False load a NN trained on 784-5-5-1
fromScratch = false;

%% Training Data
load('mnist.mat')

images = training.images;
labels = training.labels;

% Utilizing only the 1s and 7s
ones_sevens_idx = logical((labels==7) + (labels==1));
images = images(:,:,ones_sevens_idx);
labels = labels(ones_sevens_idx);

% Training Parameters
[x,y,ntrain] = size(images);
input_train = reshape(images, [x*y,ntrain]);
output_train = zeros(2,ntrain);             % Was changed from 1 to 2?
for i=1:ntrain
    if labels(i) == 1
        output_train(:,i) = [1;0];
    else
        output_train(:,i) = [0;1];
    end
end

%% Testing Data
images = test.images;
labels = test.labels;

% keep only 1's and 7's
ones_sevens_idx = logical((labels==7) + (labels==1));
images = images(:,:,ones_sevens_idx);
labels = labels(ones_sevens_idx);

[x,y,ntest] = size(images);
input_test = reshape(images, [x*y,ntest]);
output_test = zeros(2,ntest);               % Was changed from 1 to 2?
for i=1:ntest
    if labels(i) == 1
        output_test(:,i) = [1;0];
    else
        output_test(:,i) = [0;1];
    end
end

%% Training Network
if(fromScratch)
    layers = [featureInputLayer(x*y) 
             fullyConnectedLayer(5)
             reluLayer
             fullyConnectedLayer(5)
             reluLayer
             fullyConnectedLayer(2)             % Was changed from 1 to 2?
             regressionLayer];
         
    options = trainingOptions("sgdm",'Momentum', 0.95, MaxEpochs = 100,Plots='training-progress');
         
    net = trainNetwork(input_train', output_train', layers, options);

    save('mnist_784_5_5_2_new.mat')
else
    load('mnist_784_5_5_2_new.mat')
end

%% Validate the Network

output_test2 = predict(net,input_test')';
prediction_error = vecnorm( output_test - output_test2);              % Norm across each column

output_of_1s = output_test2(:,labels==1);
prediction_error_1s = prediction_error(labels==1);

output_of_7s = output_test2(:,labels==7);
prediction_error_7s = prediction_error(labels==7);

worst_1 = images(:,:,prediction_error==max(prediction_error_1s));
worst_7 = images(:,:,prediction_error==max(prediction_error_7s));

% Plotting best and worst
figure('Position',[100,100,1500,600])
subplot(1,3,1)
plot( output_of_1s(1,:), output_of_1s(2,:), 'r.', 'MarkerSize',10);
hold on;

plot( output_test2(1,prediction_error==max(prediction_error_1s)), output_test2(2,prediction_error==max(prediction_error_1s)), 'r*', 'MarkerSize',10,'LineWidth',3);

plot( output_of_7s(1,:), output_of_7s(2,:), 'b.', 'MarkerSize',10);

plot( output_test2(1,prediction_error==max(prediction_error_7s)), output_test2(2,prediction_error==max(prediction_error_7s)), 'b*', 'MarkerSize',10,'LineWidth',3);

plot( 1,0, 'kx',  'MarkerSize',10, 'LineWidth', 3);
plot( 0,1, 'ks',  'MarkerSize',10, 'LineWidth', 3);

hold off; axis([-0.1,1.1,-0.1,1.1]); axis square;
legend('Test 1s','Worst 1','Test 7s','Worst 7','True 1s','True 7s');
xlabel('Confidence 1')
ylabel('Confidence 7')
title('MNIST Test Data')


subplot(1,3,2)
imshow(worst_1(:,:,1))
title('Worst 1')

subplot(1,3,3)
% In the training data, mnist_784_5_5_2 there is a five way tie for worst
% as evident in the variable defined by 'worst_7' we will only ever look at
% the first one of each.
imshow(worst_7(:,:,1)) 
title('Worst 7')

%% Plot of all worsts
% figure('Name','Worst 1s')
% for i = 1 : size(worst_1,3)
%     subplot(1,size(worst_1,3),i)        % Gets the cells of the worst 7
%     imshow(worst_1(:,:,i)) 
% end
% sgtitle('Worst 1')
% 
figure('Name','Worst 7s')
for i = 1 : size(worst_7,3)
    subplot(1,size(worst_7,3),i)        % Gets the cells of the worst 1
    imshow(worst_7(:,:,i)) 
end
sgtitle('Worst 7')

%% Converting NN to zonotopes



Ws = {double(net.Layers(2).Weights), double(net.Layers(4).Weights), double(net.Layers(6).Weights)};
bs = {double(net.Layers(2).Bias), double(net.Layers(4).Bias), double(net.Layers(6).Bias)};
