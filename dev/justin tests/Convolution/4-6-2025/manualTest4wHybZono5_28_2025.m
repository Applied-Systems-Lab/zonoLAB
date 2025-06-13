clear
close all
clc
% %
% a = 200;
% b = 28;
% c = 28;
% tic
% convKey = @(w,i,j) sprintf('C_%d_%d_%d',w,i,j)
% convKeys = @(w,I1s,I2s,I1e,I2e) arrayfun(@(i,j) convKey(w,i,j), kron(I1s:I1e,ones(1,I2e-I2s+1)),repmat(I2s:I2e,[1,I1e-I1s+1]),UniformOutput=false)
% for i = 1:a, tmp_{i} = convKeys(i,1,1,b,c); end; out = [tmp_{:}];
% toc
% 
% tic
% for i = 1:a, tmp2_{i} = genMatIndexKeysRColumnRow(sprintf('C_%d',i),1,1,b,c); end; out2 = [tmp2_{:}];
% toc


%%
load("mnist_c5x5_20_p.mat")

%% Setting up test data point

% Defined as the average of all ones data points
weights = net1.Layers(2,1).Weights;
bias = net1.Layers(2,1).Bias;

% Collects all ones from data set puts it into a 3D matrix (28-28-xxxx)
ones_data = cat(3,XTest(:,:,XTestLabels==1), XTrain(:,:,XTrainLabels==1), XValidation(:,:,XValidLabels==1));

% Input Size
inputHeight = size(ones_data,1)
inputWidth = size(ones_data,2)

ones_data2D = reshape(ones_data,[size(ones_data,1)*size(ones_data,2) size(ones_data,3)]);

% Takes the average which will become the center point
avg_ones_data = reshape(mean(ones_data2D,2),[28 28]);


%% HybZono method, current verification lies in using the center
%% Calculating using hybzono
% in actual implmentation no need to store every cell

weights = net1.Layers(2,1).Weights;
bias = net1.Layers(2,1).Bias;

% Assumes all convolutional layers are the same if not need to apply an
% iterate
convLength = size(weights(:,:,1,1),2);
convHeight = size(weights(:,:,1,1),1);

% Creates empty matrix accounting for padding?
Empty = zeros(size(avg_ones_data,1)+2*(convHeight-1),size(avg_ones_data,2)+2*(convLength-1));

% Corresponds to Lines 68-78
%% TOPLIZIFY FUNCTION
%% WIP: Fixed to be a few lines using topletiz
stackedConvWeights = [];
stackedConvBias = [];
tic
for i = 1: length(weights)
    k = 1;
    for updown = 1: size(Empty,1)-convHeight+1
        for leftright = 1: size(Empty,2)-convLength+1
            temp = Empty;
            temp(updown:updown+(convHeight-1),leftright:leftright+(convLength-1))= rot90(weights(:,:,1,i),2);

            convShape(k,:) = reshape(temp',[size(temp,1)*size(temp,2) 1]);
            k = k + 1;
        end
    end
    stackedConvWeights = [stackedConvWeights; rot90(convShape,2)];
end
stackedConvBias = kron(reshape(bias,20,1),ones(1024,1));

toc
disp('created conv shape')
% Transforms data point to account for padding to perform full convolution
ones_data = cat(3,XTest(:,:,XTestLabels==1), XTrain(:,:,XTrainLabels==1), XValidation(:,:,XValidLabels==1));
% z = zeros(size(ones_data,1)+2*(size(weights,1)-1),size(ones_data,2)+2*(size(weights,2)-1));
z = zeros(size(ones_data,1)+2*(size(weights,1)-1),size(ones_data,2)+2*(size(weights,2)-1),size(ones_data,3));
z(size(weights,1):size(ones_data,1)+size(weights,1)-1,size(weights,2):size(ones_data,2)+size(weights,2)-1,:) = ones_data-net1.Layers(1).Mean;
ones_data = z;

% Reshapes data points into column of 784-xxxx
ones_data = reshape(ones_data,[size(ones_data,1)*size(ones_data,2) size(ones_data,3)]);

% Transforms data point to account for padding to perform full convolution
sevens_data = cat(3,XTest(:,:,XTestLabels==7), XTrain(:,:,XTrainLabels==7), XValidation(:,:,XValidLabels==7));
% z = zeros(size(ones_data,1)+2*(size(weights,1)-1),size(ones_data,2)+2*(size(weights,2)-1));
z = zeros(size(sevens_data,1)+2*(size(weights,1)-1),size(sevens_data,2)+2*(size(weights,2)-1),size(sevens_data,3));
z(size(weights,1):size(sevens_data,1)+size(weights,1)-1,size(weights,2):size(sevens_data,2)+size(weights,2)-1,:) = sevens_data-net1.Layers(1).Mean;
sevens_data = z;

% Reshapes data points into column of 784-xxxx
sevens_data = reshape(sevens_data,[size(sevens_data,1)*size(sevens_data,2) size(sevens_data,3)]);


%% ReLu weigths and Biases
Ws = {double(net1.Layers(4).Weights), double(net1.Layers(6).Weights), double(net1.Layers(8).Weights), double(net1.Layers(10).Weights)};
bs = {double(net1.Layers(4).Bias), double(net1.Layers(6).Bias), double(net1.Layers(8).Bias), double(net1.Layers(10).Bias)};

%% Average Pooling Layer
    % Creating pooling layer as one matrix
    poolLength = net1.Layers(3).PoolSize(1);
    poolHeight = net1.Layers(3).PoolSize(2);
    strideVert = net1.Layers(3).Stride(1);
    strideHorz = net1.Layers(3).Stride(2);
    
    % Number of cells averaged
    pool = 1/(poolLength*poolHeight);
    
    k = 1;
    poolMatrix = pool*ones(net1.Layers(3).PoolSize);

%% Same Convolution partitioning:
    %% Finds the partitions for a projection relative to a subset of the overall data
    snipKeys = [];
    startH = floor(convHeight/2) + 1;
    startL = floor(convLength/2) + 1;

%% Brightening Attack set up.
inputKeys = genMatIndexKeysRRowColumn('Z',1,1,inputHeight+2*(convHeight-1),inputWidth+2*(convLength-1));
sigmas = linspace(0.1,0.5,21);
data{1} = ones_data;
data{2} = sevens_data;
tic
for dataiterate = 1:2
    for i = 1: length(sigmas)
        sigma = sigmas(i);
        X1sigma{i} = interval_hull_of_points(data{dataiterate}, sigma,inputKeys);
    end
    
    % reorders the keys so its stacked per row
    reOrderKeys = genMatIndexKeysRColumnRow('Z',1,1,inputHeight+2*(convHeight-1),inputWidth+2*(convLength-1))';
    for j = 1: length(X1sigma)
        X1sigma{j} = X1sigma{j}(reOrderKeys',:,:);
    end
    %% Convolutional Implemntation change to a function
    % Probably takes in (Weights, Bias, Convolution Type)
    % Convolution type of 
        % Full
        % Same
        % Valid
    % in actual implmentation no need to store every cell

%% GENKEY FUNCTION
    outKeys = [];
    for i = 1: length(weights)
        % Create unique keys per conv layer
        outKeys{i} = genMatIndexKeysRColumnRow(sprintf('C_%d',i),1,1,inputHeight+convHeight-1,inputWidth+convLength-1)';
    end
    
    %% Use plus to merge subset
    % New  implementation using sparsity
    stackedOutKeys = vertcat(outKeys{:});

    for j = 1: length(X1sigma)
            X1conv{j} = X1sigma{j}.transform(stackedConvBias,stackedConvWeights,reOrderKeys,stackedOutKeys');
    end
   
    
%% GENKEY FUNCTION
    for i = 1: length(weights)
        % Create unique keys per conv layer
        snipKeys{i} = genMatIndexKeysRRowColumn(sprintf('C_%d',i),startH,startL, startH+inputHeight-1,startL+inputWidth-1)';
    end
    stackedSnipKeys = vertcat(snipKeys{:});

    %% Takes the projection relative to the snip keys
    for j = 1: length(X1sigma)
            X1conv{j}  = X1conv{j}(stackedSnipKeys',:,:);
    end

    %% Different Layer of the Neural Network (Average Pooling Layer)    
%% TOPLIZIFY FUNCTION
    %% Create function that will create this in one go
    for j = 1: strideVert : inputHeight
        for i = 1: strideHorz :inputWidth
            zeroMatrix = zeros(inputHeight,inputWidth);
            zeroMatrix(j:j+poolHeight-1,i:i+poolLength-1) = zeroMatrix(j:j+poolHeight-1,i:i+poolLength-1)+poolMatrix;
    
            poolingShape(k,:) = reshape(zeroMatrix',[inputHeight*inputWidth 1])';
            k = k + 1;
        end
    end
    
%% GENKEY FUNCTION
    stackedPoolingShape = kron(eye(size(weights,4)),poolingShape);
    for i = 1:length(weights)
        poolKeys{i} = genMatIndexKeysRRowColumn(sprintf('P_%d',i),1,1, 14,14);
    end
    stackedPoolKeys = [poolKeys{:}];

%% Already Implemented, ReluNN

    for j = 1: length(X1sigma)
        xHybPoolOut1 = X1conv{j}.transform([],stackedPoolingShape,stackedSnipKeys,stackedPoolKeys);

        [NN, Z1] = reluNN(xHybPoolOut1, Ws, bs, 1000);
        % 
        % [lowerBounds(j,dataiterate),upperBounds(j,dataiterate)] = bounds(Z1.Z('y_1'));
        toc
        disp('Past Bounds')
    end
end
% %%
% figure('Name','Output Range')
% shade(sigmas,lowerBounds(:,1),'r',sigmas,upperBounds(:,1),'r','FillType',[1 2;2 1]);
% hold on
% shade(sigmas,lowerBounds(:,2),'b',sigmas,upperBounds(:,2),'b','FillType',[1 2;2 1]);
% legend('','','1s','','','7s')
% grid on;
% fontsize(gcf,scale=1.5)
% xlabel('fraction of σ')
% ylabel('prediction interval')
% xlim([min(sigmas),max(sigmas)])

%% Extra functions for memzono functionality / plotting
function [Z] = interval_hull_of_points(data,threshold,keys)
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
    Z = memZono(Z,keys);
end

function [I] = plot_image_versions(Xs,N)
    I = zeros(N*36,length(Xs)*36);
    for i = 1:N
        ii = 36*(i-1)+1:36*i;
        for j = 1:length(Xs)
            jj = 36*(j-1)+1:36*j;
            X = Xs{j};
            if i == 1
                xi = ones(size(X.Gc,2),1);
            elseif i == N
                xi = -ones(size(X.Gc,2),1);
            else
                xi = 2*randi([0,1],size(X.Gc,2),1)-1;
            end
            I(ii,jj) = reshape( X.c_ + X.Gc_*xi , [36,36] );
        end
    end
end

%% Generates Keys with I rows and J columns
function keys = genMatIndexKeysRRowColumn(key,Is, Js, Ie,Je)
    keys = {};
    for j = Js:Je
        for i = Is: Ie
            keys = [keys,memZono.genKeys(sprintf("%s_%d",key,i),j)];
        end
    end
end

%% Generates Keys with I rows and J columns
function keys = genMatIndexKeysRColumnRow(key,Is, Js, Ie,Je)
    keys = {};
    for i = Js:Je
        for j = Is: Ie
            keys = [keys,memZono.genKeys(sprintf("%s_%d",key,i),j)];
        end
    end
end