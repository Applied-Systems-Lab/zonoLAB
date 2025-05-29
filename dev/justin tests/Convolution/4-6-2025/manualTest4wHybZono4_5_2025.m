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


%% Brightening Attack set up.
inputKeys = genMatIndexKeysRRowColumn('Z',1,1,inputHeight+2*(convHeight-1),inputWidth+2*(convLength-1));
sigmas = linspace(0.1,0.5,21);
data{1} = ones_data;
data{2} = sevens_data;
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
    %%
    % in actual implmentation no need to store every cell
    outKeys = [];
    for i = 1: length(weights)
        % Create unique keys per conv layer
        outKeys{i} = genMatIndexKeysRColumnRow(sprintf('C_%d',i),1,1,inputHeight+convHeight-1,inputWidth+convLength-1)';
    end
    
    %% Use plus to merge subset
    %%
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
    stackedBias = [];
    stackedWeights = [];
    stackedOutKeys = [];
    for i = 1: length(weights)
        stackedBias = [stackedBias; stackedConvBias{i}];
        stackedWeights = [stackedWeights; rot90(stackedConvWeights{i},2)];
        stackedOutKeys = [stackedOutKeys; outKeys{i}];
    end
    
    for j = 1: length(X1sigma)
            X1conv{j} = X1sigma{j}.transform(stackedBias,stackedWeights,reOrderKeys,stackedOutKeys');
    end
    
    %% Old as of 5/27/2025
    for j = 1: length(X1sigma)
        for i = 1: length(weights)
            X1convold{i,j} = X1sigma{j}.transform(stackedConvBias{i},rot90(stackedConvWeights{i},2),reOrderKeys,outKeys{i}');
        end
    end
    


    disp('past conv calc line 268')
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
    %% Finds the partitions for a projection relative to a subset of the overall data
    snipKeys = [];
    startH = floor(convHeight/2) + 1;
    startL = floor(convLength/2) + 1;
    
    for i = 1: length(weights)
        % Create unique keys per conv layer
        snipKeys{i} = genMatIndexKeysRRowColumn(sprintf('C_%d',i),startH,startL, startH+inputHeight-1,startL+inputWidth-1)';
    end
    
    %% Takes the projection relative to the snip keys
    % some reason snipKeys need to be transposed cause dim keys are row
    % vector.
    stackedSnipKeys = [];
    for i = 1: length(weights)
        stackedSnipKeys = [stackedSnipKeys; snipKeys{i}];
    end
    
    for j = 1: length(X1sigma)
            X1conv1{j}  = X1conv{j}(stackedSnipKeys',:,:);
    end


%% Old as of 5/27/2025
    % Should specifically match layerProjection(1,1)
    for j = 1: length(X1sigma)
        for i = 1: length(weights)
            X1conv1old{i,j}  = X1convold{i,j}(snipKeys{i}',:,:);
        end
    end
    
    %% 5/27 above are equal to one another maybe
    % Corresponds to Lines 94-124
    %% Average Pooling Layer
    %% Creating pooling layer as one matrix
    poolLength = net1.Layers(3).PoolSize(1);
    poolHeight = net1.Layers(3).PoolSize(2);
    strideVert = net1.Layers(3).Stride(1);
    strideHorz = net1.Layers(3).Stride(2);
    
    % Number of cells averaged
    pool = 1/(poolLength*poolHeight);
    
    k = 1;
    poolMatrix = pool*ones(net1.Layers(3).PoolSize);
    
    for j = 1: strideVert : inputHeight
        for i = 1: strideHorz :inputWidth
            zeroMatrix = zeros(inputHeight,inputWidth);
            zeroMatrix(j:j+poolHeight-1,i:i+poolLength-1) = zeroMatrix(j:j+poolHeight-1,i:i+poolLength-1)+poolMatrix;
    
            poolingShape(k,:) = reshape(zeroMatrix',[inputHeight*inputWidth 1])';
            k = k + 1;
        end
    end
    
    %% Pooling Calculation
    for i = 1: length(weights)
            % Create unique keys per conv layer
            % change pooling key definition hard coded
            poolKeys{i} = genMatIndexKeysRRowColumn(sprintf('P_%d',i),1,1, 14,14);
    end
    

    %% old as of 5/27/2025
    for j = 1: length(X1sigma)
        for i = 1: length(weights)
            xHybPoolOut1old{i,j} = X1conv1old{i,j}.transform(zeros(196,1),poolingShape,snipKeys{i},poolKeys{i});
        end
    end

    %%
    for j = 1: 1
        xConvOut1{j} = plus(xHybPoolOut1old{1,j},xHybPoolOut1old{2,j});
        for i = 3:20
           xConvOut1{j} = plus(xConvOut1{j},xHybPoolOut1old{i,j});
        end
    end
    
    disp('past combine calc')

    %% Working as of 4-6-2025
    
    % Combine each memzono individually to form giant memzono
    %%
    newPoolKeys = [];
    for i = 1:20
        newPoolKeys = [newPoolKeys poolKeys{i}];
    end
    
    disp('past pool calc')
    %%
    for j = 1: 1
        X1{j} = xConvOut1{j}(newPoolKeys,:,:);
    end
    % Corresponds to Lines 126-170
%%
    stackedPoolingShape = kron(eye(size(weights,4)),poolingShape);
    stackedPoolKeys = [];
    for i = 1: length(weights)
        stackedPoolKeys = [stackedPoolKeys; poolKeys{i}'];
    end
%%
    for j = 1: length(X1sigma)
            xHybPoolOut1{j} = X1conv1{j}.transform(zeros(size(stackedPoolingShape,1),1),stackedPoolingShape,stackedSnipKeys,stackedPoolKeys');
    end


    % Taking ReLU with ReLUNN 
    % First Relu Layer (xxxx*50)
    % W * x + B
    Ws = {double(net1.Layers(4).Weights), double(net1.Layers(6).Weights), double(net1.Layers(8).Weights), double(net1.Layers(10).Weights)};
    bs = {double(net1.Layers(4).Bias), double(net1.Layers(6).Bias), double(net1.Layers(8).Bias), double(net1.Layers(10).Bias)};
    
    %%
    for j = 1: length(X1sigma)
        [NN1, Z1] = reluNN(xHybPoolOut1{j}, Ws, bs, 1000);
        
        % Current issue with bounds
        temp = Z1.Z('y_1');
        temp.b = double(temp.b);
        [lb1{j,dataiterate},ub1{j,dataiterate}] = temp.bounds;
        disp('past bounds')
    end
end
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