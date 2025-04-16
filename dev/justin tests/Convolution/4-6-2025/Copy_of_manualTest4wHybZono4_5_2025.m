clear
close all
clc

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

%% Step by step calculation through NN
% Stores expected values for the any projection
testPoint = avg_ones_data;

newNet = dlnetwork;
layerProjection = [];

layers = [
    imageInputLayer([28 28 1],Mean = net1.Layers(1).Mean)]
newNet = addLayers(newNet,layers);
newNet = initialize(newNet)
newNet.predict(testPoint)

for i = 1: length(net1.Layers)-1
    layers = net1.Layers(i+1)
    newNet = addLayers(newNet,layers)
    newNet = connectLayers(newNet, newNet.Layers(i).Name, newNet.Layers(i+1).Name)
    newNet = initialize(newNet)
    layerProjection{i} = newNet.predict(testPoint)
end

%% Propagating through the network step by step
%% Convolution Layer Projection

convLength = size(weights(:,:,1,i),2);
convHeight = size(weights(:,:,1,i),1);
Empty = zeros(size(testPoint,1)+2*(convHeight-1),size(testPoint,2)+2*(convLength-1));

%% Old manual Calculation Process
testPoint = testPoint-net1.Layers(1).Mean;

Ashape = Empty;
Ashape(convHeight:size(Ashape,1)-convHeight+1,convLength:size(Ashape,2)-convLength+1) = testPoint;
Ashape = reshape(Ashape',[size(Ashape,1)*size(Ashape,2) 1]);
for i = 1: length(weights)
    % Generates a single matrix to calculate the convolution
    k = 1;
    for updown = 1: size(Empty,1)-convHeight+1
        for leftright = 1: size(Empty,2)-convLength+1
            temp = Empty;
            temp(updown:updown+(convHeight-1),leftright:leftright+(convLength-1))= rot90(weights(:,:,1,i),2);

            convShape(k,:) = reshape(temp',[size(temp,1)*size(temp,2) 1]);
            k = k + 1;
        end
    end
    oldManualconvShape{i,1} = convShape;
    oldManualbias(i,1) = bias(:,:,i);

    % Full convolution
    convSquare = reshape(rot90(convShape,2)*Ashape+bias(:,:,i),[size(testPoint,2)+convLength-1 size(testPoint,1)+convHeight-1])';
    oldManualconvSquare{i,1} = convSquare;      % Should match with full conv below [convSquare]
    % Same Convolution partitioning:
    startH = floor(convHeight/2) + 1;
    startL = floor(convLength/2) + 1;

    forward{i} = convSquare(startH:startH+size(testPoint,1)-1,startL:startL+size(testPoint,2)-1);
end

% Reshaping Input for Pooling
x=[];
for i = 1:length(net1.Layers(2).Weights)
    x = [x ;reshape(forward{i},[1 size(forward{i},1)*size(forward{i},2)])];
end

% Average Pooling Layer calculation
% Assigning pooling parameters
poolLength = net1.Layers(3).PoolSize(1);
poolHeight = net1.Layers(3).PoolSize(2);
strideVert = net1.Layers(3).Stride(1);
strideHorz = net1.Layers(3).Stride(2);

% Number of cells averaged
pool = 1/(poolLength*poolHeight);

k = 1;
poolMatrix = pool*ones(net1.Layers(3).PoolSize);

% Generates a single matrix to calculate pooling
for j = 1: strideVert : size(forward{1},1)
    for i = 1: strideHorz :size(forward{1},2)
        zeroMatrix = zeros(size(forward{1}));
        zeroMatrix(j:j+poolHeight-1,i:i+poolLength-1) = zeroMatrix(j:j+poolHeight-1,i:i+poolLength-1)+poolMatrix;

        poolingShape(k,:) = reshape(zeroMatrix',[size(forward{1},1)*size(forward{1},2) 1])';
        k = k + 1;
    end
end

xPoolOut = poolingShape*x';    % calcs fine

% RELU Propagation
%% Relu layer 4-5
weights = net1.Layers(4,1).Weights;
xPoolOut = reshape(xPoolOut,[size(weights,2) 1]);
bias = net1.Layers(4,1).Bias;

% linear mapping
% Validating layerProjection(1,3)
forwardRelu1 = weights*xPoolOut+bias;

% relu activation function
% Validating layerProjection(1,4)
activeRelu1 = max(0,forwardRelu1);

%% Relu layer 6-7
weights = net1.Layers(6,1).Weights;
bias = net1.Layers(6,1).Bias;

% linear mapping
% Validating layerProjection(1,5)
forwardRelu2 = weights*activeRelu1 + bias;

% relu activation function
% Validating layerProjection(1,6)
activeRelu2 = max(0,forwardRelu2);

%% Relu layer 8-9
weights = net1.Layers(8,1).Weights;
bias = net1.Layers(8,1).Bias;

% linear mapping
% Validating layerProjection(1,7)
forwardRelu3 = weights*activeRelu2 + bias;

% relu activation function
% Validating layerProjection(1,8)
activeRelu3 = max(0,forwardRelu3);

%% Relu layer 10
weights = net1.Layers(10,1).Weights;
bias = net1.Layers(10,1).Bias;

% linear mapping
% Validating layerProjection(1,9)
lastLayer = weights*activeRelu3 + bias


%% HybZono method, current verification lies in using the center
%% Calculating using hybzono
% in actual implmentation no need to store every cell

weights = net1.Layers(2,1).Weights;
bias = net1.Layers(2,1).Bias;

% Corresponds to Lines 68-78
stackedConvWeights = [];
stackedConvBias = [];
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
    stackedConvWeights{i,1}  = convShape;
    stackedConvBias{i,1}    = bias(:,:,i)*ones(1024,1);
end
% % % Compare: oldManualconvShape to stackedConvWeights
% % for i = 1: 200
% %     if(norm(oldManualconvShape{i} - stackedConvWeights{i}))
% %         disp(norm(oldManualconvShape{i} - stackedConvWeights{i}))
% %     end
% % end
% % 
% % % Compare: oldManualbias to stackedConvBias
% % for i = 1: 200
% %     if(norm(oldManualbias(i) - stackedConvBias{i}(1,1)))
% %         disp(norm(oldManualbias(i) - stackedConvBias{i}(1,1)))
% %     end
% % end

% Transforms data point to account for padding to perform full convolution
ones_data = cat(3,XTest(:,:,XTestLabels==1), XTrain(:,:,XTrainLabels==1), XValidation(:,:,XValidLabels==1));
z = zeros(size(ones_data,1)+2*(size(weights,1)-1),size(ones_data,2)+2*(size(weights,2)-1),size(ones_data,3));
z(size(weights,1):size(ones_data,1)+size(weights,1)-1,size(weights,2):size(ones_data,2)+size(weights,2)-1,:) = ones_data-net1.Layers(1).Mean;
ones_data = z;

% Reshapes data points into column of 784-xxxx
ones_data = reshape(ones_data,[size(ones_data,1)*size(ones_data,2) size(ones_data,3)]);

% Transforms data point to account for padding to perform full convolution
sevens_data = cat(3,XTest(:,:,XTestLabels==7), XTrain(:,:,XTrainLabels==7), XValidation(:,:,XValidLabels==7));
z = zeros(size(sevens_data,1)+2*(size(weights,1)-1),size(sevens_data,2)+2*(size(weights,2)-1),size(sevens_data,3));
z(size(weights,1):size(sevens_data,1)+size(weights,1)-1,size(weights,2):size(sevens_data,2)+size(weights,2)-1,:) = sevens_data-net1.Layers(1).Mean;
sevens_data = z;

% Reshapes data points into column of 784-xxxx
sevens_data = reshape(sevens_data,[size(sevens_data,1)*size(sevens_data,2) size(sevens_data,3)]);

% data variable for later use
data{1} = ones_data;
data{2} = sevens_data;

%% %% Key Generation for ordering
%% Initial Set Up
% creates key order for input vector (counting rows then columns)
inputKeys = genMatIndexKeysRRowColumn('Z',1,1,inputHeight+2*(convHeight-1),inputWidth+2*(convLength-1));

% reorders the keys for Xsigma to get it into the form of (counting columns then rows)
reOrderKeys = genMatIndexKeysRColumnRow('Z',1,1,inputHeight+2*(convHeight-1),inputWidth+2*(convLength-1))';

%% Conv Key Set Up
% Key set up for the convolutional kernels (counting columns then rows)
outKeys = [];
for i = 1: length(weights)
    % Create unique keys per conv layer
    outKeys{i} = genMatIndexKeysRColumnRow(sprintf('C_%d',i),1,1,inputHeight+convHeight-1,inputWidth+convLength-1)';
end

%% Finds the partitions for a projection relative to a subset of the overall data
snipKeys = [];
startH = floor(convHeight/2) + 1;
startL = floor(convLength/2) + 1;

% Keys refering to a subset of the full matrix (counting rows then columns)
for i = 1: length(weights)
    % Create unique keys per conv layer
    snipKeys{i} = genMatIndexKeysRRowColumn(sprintf('C_%d',i),startH,startL, startH+inputHeight-1,startL+inputWidth-1)';
end

%% Average Pooling Layer
% Creating pooling layer as one matrix

% Definition of pooling parameters
poolLength = net1.Layers(3).PoolSize(1);
poolHeight = net1.Layers(3).PoolSize(2);
strideVert = net1.Layers(3).Stride(1);
strideHorz = net1.Layers(3).Stride(2);

% Number of cells averaged
pool = 1/(poolLength*poolHeight);

k = 1;
poolMatrix = pool*ones(net1.Layers(3).PoolSize);
% Single matrix for pooling (counting rows then columns)
for i = 1: length(weights)
        % Create unique keys per conv layer
        % change pooling key definition hard coded
        poolKeys{i} = genMatIndexKeysRRowColumn(sprintf('P_%d',i),1,1, 14,14);
end

% Stacks pool keys for combining layer to produce single hybzono out
newPoolKeys = [];
for i = 1:20
    newPoolKeys = [newPoolKeys poolKeys{i}];
end

% Taking ReLU with ReLUNN 
% Ws * x + bs
Ws = {double(net1.Layers(4).Weights), double(net1.Layers(6).Weights), double(net1.Layers(8).Weights), double(net1.Layers(10).Weights)};
bs = {double(net1.Layers(4).Bias), double(net1.Layers(6).Bias), double(net1.Layers(8).Bias), double(net1.Layers(10).Bias)};

%% Brightening Attack Calculation
sigmas = linspace(0.1,0.5,21);
tic
for dataiterate = 1:2
            tic
    for i = 1: length(sigmas)
        sigma = sigmas(i);
        Xsigma{i} = interval_hull_of_points(data{dataiterate}, sigma,inputKeys);
    end
    
    for j = 1: length(Xsigma)
        Xsigma{j} = Xsigma{j}(reOrderKeys',:,:);
    end

    % Corresponds to Lines 85
    % Should match up with [convSquare] above
    
    % for validation
    
    % input check
    % reshape(table2array(X1sigma{1,1}.c),36,36)-reshape(Ashape,[36 36]) % this compares if the input is the same as the oldmanual one
    % as of 4-6-2025 difference is negligible 
    
    % weights check
    % norm(rot90(oldManualconvShape{1},2) - rot90(stackedConvWeights{1},2))
    % as of 4-6-2025 no difference
    
    % bias check
    % norm(bias(:,:,1)*ones(1024,1)-stackedConvBias{1})
    % as of 4-6-2025 no difference
    
    % calculation of convSquare from old
    % convSquare = reshape(rot90(oldManualconvShape{i},2)*Ashape+bias(:,:,i),[size(testPoint,2)+convLength-1 size(testPoint,1)+convHeight-1])';
    for j = 1: length(Xsigma)
        for i = 1: length(weights)
            Xconv{i,j} = Xsigma{j}.transform(stackedConvBias{i},rot90(stackedConvWeights{i},2),reOrderKeys,outKeys{i}');
        end
    end
            toc
    disp('past conv calc line 268')
            tic
    % 
    % for i = 1: length(weights)
    %     reOrderKeysConv = genMatIndexKeysRRowColumn(sprintf('C_%d',i),1,1,inputHeight+convHeight-1,inputWidth+convLength-1)';
    %     X1conv1{i,1} = X1conv1{i,1}(reOrderKeysConv',:,:);
    % end
    % resolved utilizing wrong input keys, swapped to reorder and obtained the
    % proper output for case of 1 input image. Not tested on full data set.
    
    % This is the hybzono representation of the oldmanualconvsquare
    % conv1matrix = reshape(table2array(X1conv1{1,1}.c),[32 32])';
    
    % Corresponds to Lines 87-91
    %% Same Convolution partitioning:
    %% Takes the projection relative to the snip keys
    % Should specifically match layerProjection(1,1)
    for j = 1: length(Xsigma)
        for i = 1: length(weights)
            XconvOrdered{i,j}  = Xconv{i,j}(snipKeys{i}',:,:);
        end
    end
    
    % Corresponds to Lines 94-124
    for j = 1: strideVert : inputHeight
        for i = 1: strideHorz :inputWidth
            zeroMatrix = zeros(inputHeight,inputWidth);
            zeroMatrix(j:j+poolHeight-1,i:i+poolLength-1) = zeroMatrix(j:j+poolHeight-1,i:i+poolLength-1)+poolMatrix;
    
            poolingShape(k,:) = reshape(zeroMatrix',[inputHeight*inputWidth 1])';
            k = k + 1;
        end
    end
    
    for j = 1: length(Xsigma)
        for i = 1: length(weights)
            xHybPoolOut1{i,j} = XconvOrdered{i,j}.transform(zeros(196,1),poolingShape,snipKeys{i},poolKeys{i});
        end
    end
    
    %% Working as of 4-6-2025
            toc
    disp('past pool calc')
            tic
    %%
    % Combine each memzono individually to form giant memzono
    for j = 1: length(Xsigma)
        xConvOut1{j} = plus(xHybPoolOut1{1,j},xHybPoolOut1{2,j});
        for i = 3:20
           xConvOut1{j} = plus(xConvOut1{j},xHybPoolOut1{i,j});
        end
    end
            toc
    disp('past combine calc')
            tic
    % Corresponds to Lines 126-170
    %%
    for j = 1: length(Xsigma)
        X1{j} = xConvOut1{j}(newPoolKeys,:,:);
    end
    
    %%
    for j = 1: length(Xsigma)
        [NN1, Z1] = reluNN(X1{j}, Ws, bs, 1000);
        
        %%
        % Current issue with bounds
        temp = Z1.Z('y_1');
        temp.b = double(temp.b);
        [lb1{j,dataiterate},ub1{j,dataiterate}] = temp.bounds;
    end
            toc
    disp('past bounds')
end
toc
%%
lowerBounds = cell2mat(lb1);
upperBounds = cell2mat(ub1);
figure('Name','Output Range')
shade(sigmas,lowerBounds(:,1),'r',sigmas,upperBounds(:,1),'r','FillType',[1 2;2 1]);
hold on
shade(sigmas,lowerBounds(:,2),'b',sigmas,upperBounds(:,2),'b','FillType',[1 2;2 1]);
legend('','','1s','','','7s')
grid on;
fontsize(gcf,scale=1.5)
xlabel('fraction of σ')
ylabel('prediction interval')
xlim([min(sigmas),max(sigmas)])

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

%% Generates Keys with I rows per J columns
function keys = genMatIndexKeysRRowColumn(key,Is, Js, Ie,Je)
    keys = {};
    for j = Js:Je
        for i = Is: Ie
        keys = [keys,memZono.genKeys(sprintf("%s_%d",key,i),j)];
        end
    end
end

%% Generates Keys with J columns per I rows
function keys = genMatIndexKeysRColumnRow(key,Is, Js, Ie,Je)
    keys = {};
    for i = Js:Je
        for j = Is: Ie
        keys = [keys,memZono.genKeys(sprintf("%s_%d",key,i),j)];
        end
    end
end