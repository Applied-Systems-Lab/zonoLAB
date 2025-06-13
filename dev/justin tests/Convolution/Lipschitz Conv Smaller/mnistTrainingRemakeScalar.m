% Programmer Name: Justin Chen
% The following code is a updated version of the 2022-11-29 version of the
% mnist command file form the L4DC papers.

% List of Changes
% Changed the output_bounds to reference the updated bounds function found
% in zonoLab
% Changed the interval_hull_of_points to the current variable convention
% Added the reluNN to the @memZono folder so there is a memory aspect?

clear; close all;
clc;

% True trains the mnist data
% False load a NN trained on 784-5-5-1
fromScratch = false;

%% Training Data
load('5x5mnist.mat')

% Utilizing only the 1s and 7s
ones_zeros_idx = logical((training.labels==0) + (training.labels==1));
training.images = training.images(:,:,ones_zeros_idx);
training.labels = training.labels(ones_zeros_idx);

% Training Parameters
[x,y,ntrain] = size(training.images);
input_train = reshape(training.images, [x*y,ntrain]);
output_train = zeros(1,ntrain);             % Was changed from 1 to 2?
for i=1:ntrain
    if training.labels(i) == 1
        output_train(i) = 1;
    else
        output_train(i) = -1;
    end
end


%% Testing Data
% keep only 1's and 7's
ones_zeros_idx = logical((test.labels==0) + (test.labels==1));
test.images = test.images(:,:,ones_zeros_idx);
test.labels =  test.labels(ones_zeros_idx);

[x,y,ntest] = size(test.images);
input_test = reshape(test.images, [x*y,ntest]);
output_test = zeros(1,ntest);               % Was changed from 1 to 2?
for i=1:ntest
    if test.labels(i) == 1
        output_test(i) = 1;
    else
        output_test(i) = -1;
    end
end

%% Training Network
if(fromScratch)
    layers = [featureInputLayer(x*y) 
             fullyConnectedLayer(5)
             reluLayer
             fullyConnectedLayer(5)
             reluLayer
             fullyConnectedLayer(1)             % Was changed from 1 to 2?
             regressionLayer];
         
    options = trainingOptions("sgdm",'Momentum', 0.95, MaxEpochs = 100,Plots='training-progress');
         
    net = trainNetwork(input_train', output_train', layers, options);

    save('mnist_25_5_5_1_new.mat', 'net')
else
    % Newly trained use 'mnist_784_5_5_1_new.mat'
    % L4DC Paper trained on 'mnist_784_5_5_1_all.mat'
    load('mnist_25_5_5_1_new.mat','net')
end

%% Validate the Network

output_test2 = predict(net,input_test')';
prediction_error = abs( output_test - output_test2);              % Norm across each column

classification = sign(output_test2);
accuracy = sum(classification==output_test)/length(classification)*100;

output_of_1s = output_test2(test.labels==1);
prediction_error_1s = prediction_error(test.labels==1);

output_of_0s = output_test2(test.labels==0);
prediction_error_0s = prediction_error(test.labels==0);

worst_1 = test.images(:,:,prediction_error==max(prediction_error_1s));
worst_0 = test.images(:,:,prediction_error==max(prediction_error_0s));

% Plotting best and worst
figure('Name','HistoGram','Position',[100,100,1500,600])
subplot(1,3,1)
histogram( output_of_1s, 'Normalization','probability' );
hold on;
histogram( output_of_0s, 'Normalization','probability' );
plot( 1,0, 'kx',  'MarkerSize',10, 'LineWidth', 3);
plot( -1,0, 'ks',  'MarkerSize',10, 'LineWidth', 3);
hold off;
legend('Test 1s','Test 0s','True 1s','True 0s');

xlabel('Number')
ylabel('Confidence')
title('MNIST Test Data')

subplot(1,3,2)
imshow(worst_1(:,:,1))
title('Worst 1')

subplot(1,3,3)
% In the training data, mnist_784_5_5_2 there is a five way tie for worst
% as evident in the variable defined by 'worst_7' we will only ever look at
% the first one of each.
imshow(worst_0(:,:,1)) 
title('Worst 0')

% Plot of all worsts
% If there are multiple worsts
figure('Name','Worst 1s')
for i = 1 : size(worst_1,3)
    subplot(1,size(worst_1,3),i)        % Gets the cells of the worst 7
    imshow(worst_1(:,:,i)) 
end
sgtitle('Worst 1')

figure('Name','Worst 0s')
for i = 1 : size(worst_0,3)
    subplot(1,size(worst_0,3),i)        % Gets the cells of the worst 1
    imshow(worst_0(:,:,i)) 
end
sgtitle('Worst 0')

%% Converting NN to zonotopes

% From the original paper, interval_hull_of_points constructs a zono obj so
% it does not need to be changed.

%% Define Input Sets
X = interval_hull_of_points([input_train, input_test], 'all');

ones_data = [ input_train(:,training.labels==1) , input_test(:,test.labels==1) ];
zeros_data = [ input_train(:,training.labels==0) , input_test(:,test.labels==0) ];
X_ones = interval_hull_of_points( ones_data, 'all' ); % interval Hull of all ones images
X_zeros = interval_hull_of_points( zeros_data, 'all' ); % interval Hull of all sevens images

obvious_ones_data = keep_closest(ones_data,250);
obvious_zeros_data = keep_closest(zeros_data,250);
X_obvious_ones = interval_hull_of_points( obvious_ones_data, 1 );
X_obvious_zeros = interval_hull_of_points( obvious_zeros_data, 1 );
X_obvious_ones_good = interval_hull_of_points( obvious_ones_data, 0.5 ); %0.1
X_obvious_zeros_good = interval_hull_of_points( obvious_zeros_data, 0.5 ); %0.1

%% Average of set data
figure('Name','Overlaided inputs')

% Collage of all inputs
subplot(2,3,1)
imshow(reshape( table2array(X.c) , [5,5] ))
title('all input set')

% Collage of all 1s
subplot(2,3,2)
imshow(reshape( table2array(X_ones.c) , [5,5] ))
title('all ones')
subplot(2,3,3)

% Collage of the obvious 1s
imshow(reshape( table2array(X_obvious_ones.c) , [5,5] ))
title('obvious ones')

% Collage of all 7s
subplot(2,3,4)
imshow(reshape( table2array(X_zeros.c) , [5,5] ))
title('all zeros')


% Collage of the obvious 7s
subplot(2,3,5)
imshow(reshape( table2array(X_obvious_zeros.c) , [5,5] ))
title('obvious zeros')

%% Construct Hybrid Zonotopes
% Ws and bs are hardcoded and might need to change if the layers are
% manipulated
Ws = {double(net.Layers(2).Weights), double(net.Layers(4).Weights), double(net.Layers(6).Weights)};
bs = {double(net.Layers(2).Bias), double(net.Layers(4).Bias), double(net.Layers(6).Bias)};

fprintf('All 1s and 0s:\n')
Z_range = output_bound(X,Ws,bs);
fprintf('Range is [ %1.4f , %1.4f ]\n\n', Z_range);

fprintf('All 1s:\n')
Z_ones_range = output_bound(X_ones,Ws,bs);
fprintf('Range is [ %1.4f , %1.4f ]\n\n', Z_ones_range);

fprintf('All 0s:\n')
Z_zeros_range = output_bound(X_zeros,Ws,bs);
fprintf('Range is [ %1.4f , %1.4f ]\n\n', Z_zeros_range);

fprintf('All Obvious 1s:\n')
Z_obvious_ones_range = output_bound(X_obvious_ones,Ws,bs);
fprintf('Range is [ %1.4f , %1.4f ]\n\n', Z_obvious_ones_range);

fprintf('All Obvious 0s:\n')
Z_obvious_zeros_range = output_bound(X_obvious_zeros,Ws,bs);
fprintf('Range is [ %1.4f , %1.4f ]\n\n', Z_obvious_zeros_range);
%%
fprintf('All Obvious (Good) 1s:\n')
Z_obvious_ones_good_range = output_bound(X_obvious_ones_good,Ws,bs);
fprintf('Range is [ %1.4f , %1.4f ]\n\n', Z_obvious_ones_good_range);

fprintf('All Obvious (Good) 0s:\n')
Z_obvious_zeros_good_range = output_bound(X_obvious_zeros_good,Ws,bs);
fprintf('Range is [ %1.4f , %1.4f ]\n\n', Z_obvious_zeros_good_range);

%% 
subplot(2,3,6)
plot( [1,1],Z_ones_range, 'r', 'LineWidth', 3); hold on;
plot( [1,1],Z_obvious_ones_range, 'b', 'LineWidth', 2);
plot( [7,7],Z_zeros_range, 'm', 'LineWidth', 3);
plot( [7,7],Z_obvious_zeros_range, 'c', 'LineWidth', 2);
plot( [4,4],Z_range, 'k', 'LineWidth', 6); hold off;
legend('1s','Obvious 1s','0s','Obvious 0s','1s and 0s')
xticks([1,4,7])
xticklabels({'1s','1s and 0s','0s'})
xlim([0,8])
ylabel('Deviation from normal')

%% Fuzziness from generators

figure('Name','Adding in uncertainty?')
subplot(2,3,1)
imshow(reshape( table2array(X.c + table2array(X.Gc)*ones(size(X.Gc,2),1)) , [5,5] ))
title('all input set')
subplot(2,3,2)
imshow(reshape( table2array(X_ones.c + table2array(X_ones.Gc)*ones(size(X_ones.Gc,2),1)) , [5,5] ))
title('all ones')
subplot(2,3,3)
imshow(reshape( table2array(X_obvious_ones.c + table2array(X_obvious_ones.Gc)*ones(size(X_obvious_ones.Gc,2),1)) , [5,5] ))
title('obvious ones')
subplot(2,3,4)
imshow(reshape( table2array(X_zeros.c + table2array(X_zeros.Gc)*ones(size(X_zeros.Gc,2),1)) , [5,5] ))
title('all sevens')
subplot(2,3,5)
imshow(reshape( table2array(X_obvious_zeros.c + table2array(X_obvious_zeros.Gc)*ones(size(X_obvious_zeros.Gc,2),1)) , [5,5] ))
title('obvious sevens')

%% Final Plots

figure('Name','Final Plot')
boxchart(1+0*output_of_1s, output_of_1s, 'BoxEdgeColor','r','MarkerColor','r');
hold on;
boxchart(0+0*output_of_0s, output_of_0s, 'BoxEdgeColor','b','MarkerColor','b');
plot([0.5,7.5],[0,0],'k--');
hold off;

xticks([1,7]);
xticklabels({'1s (test)','0s (test)'})
yticks([-1,0,1]);

inset_axes = axes('Position', [0.34 0.14 0.35 0.35]);
imshow(worst_1)
pbaspect([1 1 1])

inset_axes = axes('Position', [0.34 0.55 0.35 0.35]);
imshow(worst_0)
pbaspect([1 1 1])

%%
X1sigma = {};
X7sigma = {};
sigmas = linspace(0.1,0.5,21);
X1sigma_range = zeros(2,length(sigmas));
X7sigma_range = zeros(2,length(sigmas));
for i = 1:length(sigmas)
    sigma = sigmas(i)
    X1sigma{i} = interval_hull_of_points( ones_data, sigma );
    X7sigma{i} = interval_hull_of_points( zeros_data, sigma );
    
    X1sigma_range(:,i) = output_bound(X1sigma{i},Ws,bs);
    [X1sigma_NN{i}, X1sigma_Z{i}] = reluNN(X1sigma{i}, Ws, bs, 1000);
    X7sigma_range(:,i) = output_bound(X7sigma{i},Ws,bs);
    [X0sigma_NN{i}, X0sigma_Z{i}] = reluNN(X7sigma{i}, Ws, bs, 1000);
end

%%
figure('Name','Output Range')
shade(sigmas,X1sigma_range(1,:),'r',sigmas,X1sigma_range(2,:),'r','FillType',[1 2;2 1]);
hold on;
shade(sigmas,X7sigma_range(1,:),'b',sigmas,X7sigma_range(2,:),'b','FillType',[1 2;2 1]);
legend('','','1s','','','0s')
grid on;
fontsize(gcf,scale=1.5)
xlabel('fraction of σ')
ylabel('prediction interval')
xlim([min(sigmas),max(sigmas)])
ylim([-1.4,1.4])
yticks([-1,0,1])

%%
figure()
subplot(2,1,1)
imshow(plot_image_versions(X1sigma,5))
subplot(2,1,2)
imshow(plot_image_versions(X7sigma,5))
xlabel('σ')
ylabel('ξ=1, ξ=-1, ξ∈\{-1,1\}^{784}')

%%
figure()
xi = 2*randi([0,1],size(X.Gc,2),1)-1;
for i=1:5
subplot(1,5,i)
Xsi = X7sigma{abs(sigmas-i/10)<0.0001};
imshow( reshape( table2array(Xsi.c + table2array(Xsi.Gc)*xi) , [5,5] ) )
end

% computes the interval hull of a set of points 
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

function [core] = keep_closest(data, num)
    nominal = mean(data,2);
    error_norm = vecnorm(data-nominal,1); % evaluate how far away each image is from the nominal
    [esorted,I] = sort(error_norm);
    core = data(:,I);
    core = core(:,1:num);
end

function [range] = output_bound(X,Ws,bs)
    tic
    [NN, Z] = reluNN(X, Ws, bs, 1000);
    toc
    % For the specific example, dim order does not matter but could
    % possibly be an issue
    [range(1),range(2)] = bounds(Z.Z(Z.dimKeys));
    toc
end

function [I] = plot_image_versions(Xs,N)
    I = zeros(N*5,length(Xs)*5);
    for i = 1:N
        ii = 5*(i-1)+1:5*i;
        for j = 1:length(Xs)
            jj = 5*(j-1)+1:5*j;
            X = Xs{j};
            if i == 1
                xi = ones(size(X.Gc,2),1);
            elseif i == N
                xi = -ones(size(X.Gc,2),1);
            else
                xi = 2*randi([0,1],size(X.Gc,2),1)-1;
            end
            I(ii,jj) = reshape( table2array(X.c + table2array(X.Gc)*xi) , [5,5] );
        end
    end
end