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
% newNet = dlnetwork;
% layerProjection = [];
% 
% layers = [
%     imageInputLayer([28 28 1],Mean = net1.Layers(1).Mean)]
% newNet = addLayers(newNet,layers);
% newNet = initialize(newNet)
% newNet.predict(m7)
% 
% for i = 1: length(net1.Layers)-1
%     layers = net1.Layers(i+1)
%     newNet = addLayers(newNet,layers)
%     newNet = connectLayers(newNet, newNet.Layers(i).Name, newNet.Layers(i+1).Name)
%     newNet = initialize(newNet)
%     layerProjection{i} = newNet.predict(m7)
%     % plot(newNet)
% end
% 
% m7 = reshape(XTest(:,:,1,1),[28 28])-net1.Layers(1).Mean;
% m1 = reshape(XTest(:,:,1,2),[28 28])-net1.Layers(1).Mean; %

%%
% figure
% imshow(m7)
% figure
% imshow(m1)

%% Propagating through the network step by step
start=1;
% for j = 1:2163
testPoint = reshape(XTest(:,:,1,start),[28 28])-net1.Layers(1).Mean;
%% For the 7
weights = net1.Layers(2,1).Weights;
bias = net1.Layers(2,1).Bias;

%% Defining input zonotope
sigmas = linspace(0.1,0.5,21);
ones_data = cat(3,XTest(:,:,XTestLabels==1), XTrain(:,:,XTrainLabels==1), XValidation(:,:,XValidLabels==1))-net1.Layers(1).Mean;

z = zeros(size(testPoint,1)+2*(size(weights,1)-1),size(testPoint,2)+2*(size(weights,2)-1),size(ones_data,3));
z(size(weights,1):size(testPoint,1)+size(weights,1)-1,size(weights,2):size(testPoint,2)+size(weights,2)-1,:) = ones_data;
ones_data = z;
%%
ones_data = reshape(ones_data,[size(ones_data,1)*size(ones_data,2) size(ones_data,3)]);
%%
sevens_data = cat(3,XTest(:,:,XTestLabels==7), XTrain(:,:,XTrainLabels==7), XValidation(:,:,XValidLabels==7))-net1.Layers(1).Mean;
z = zeros(size(testPoint,1)+2*(size(weights,1)-1),size(testPoint,2)+2*(size(weights,2)-1),size(sevens_data,3));
z(size(weights,1):size(testPoint,1)+size(weights,1)-1,size(weights,2):size(testPoint,2)+size(weights,2)-1,:) = sevens_data;
sevens_data = z;
%%
sevens_data = reshape(sevens_data,[size(sevens_data,1)*size(sevens_data,2) size(sevens_data,3)]);

%%
for i = 1: length(sigmas)
    sigma = sigmas(i);
    X1sigma{i} = interval_hull_of_points(ones_data, sigma);
    X7sigma{i} = interval_hull_of_points(sevens_data, sigma);
end

%% Plot sample sigma
figure()
subplot(2,1,1)
imshow(plot_image_versions(X1sigma,5))
subplot(2,1,2)
imshow(plot_image_versions(X7sigma,5))
xlabel('σ')
ylabel('ξ=1, ξ=-1, ξ∈\{-1,1\}^{784}')

% figure
% hold on

convLength = size(weights(:,:,1,i),2);
convHeight = size(weights(:,:,1,i),1);
%% Old manual Calculation
Empty = zeros(size(testPoint,1)+2*(convHeight-1),size(testPoint,2)+2*(convLength-1));
% Ashape = Empty;
% Ashape(convHeight:size(Ashape,1)-convHeight+1,convLength:size(Ashape,2)-convLength+1) = testPoint;
% Ashape = reshape(Ashape',[size(Ashape,1)*size(Ashape,2) 1]);

stackedConvWeights = [];
stackedConvBias = [];
%% Calculating using hybzono
% in actual implmentation no need to store every cell
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
        stackedConvWeights{i}  = rot90(convShape,2);
        stackedConvBias{i}    = bias(:,:,i)*ones(1024,1);
    end

%%
% in actual implmentation no need to store every cell
    outKeys = [];
    for i = 1: length(weights)
        % Create unique keys per conv layer
        outKeys{i} = genMatIndexKeys(sprintf('C_%d',i),1,1,size(testPoint,1)+convHeight-1,size(testPoint,2)+convLength-1)';
    end

        %%
        % Full convolution calculation, output of conv layer
        % Must iterate through all sigma (21) and all weights (200)
        % convSquare{i} = X1sigma{sigIterate}.transform(bias(:,:,i)*ones(1024,1),rot90(convShape,2),memZono.genKeys('Z',1:1296),outKeys);
        
        % Use projection to obtain subset of conv
        % The convolution operations, 'same' and 'valid' are subset processes
        % of the full convolution. Applying full convolution and then
        % partitioning the proper section of the matricies will allow for
        % calculation of convolution methods other than full.
        % 'same' : the subset is specifically,
        % (size of convolution)/2  + 1 for starting location
        % then the size of the resulting conv is the size of A.
    
        % 'valid': the subset is specifically,
        % size of convolution for starting location
        % size of resulting conv is 
        %       (length of input - length of conv + 1)
        %       (width of input - width of conv + 1)

        
        % Use plus to merge subset


        % convSquare = reshape(rot90(convShape,2)*Ashape+bias(:,:,i),[size(testPoint,2)+convLength-1 size(testPoint,1)+convHeight-1])';

        %%
        % If values look weird, stackedConvWeights might be generated
        % incorrectly
    for i = 1: length(weights)
        X1conv1{i,1} = X1sigma{1}.transform(stackedConvBias{i},stackedConvWeights{i},memZono.genKeys('Z',1:1296),outKeys{i}');
    end

%% Same Convolution partitioning:
%% Finds the partitions for a projection relative to a subset of the overall data
    snipKeys = [];
    startH = floor(convHeight/2) + 1;
    startL = floor(convLength/2) + 1;

    for i = 1: length(weights)
        % Create unique keys per conv layer
        snipKeys{i} = genMatIndexKeys(sprintf('C_%d',i),startH,startL, startH+size(testPoint,1)-1,startL+size(testPoint,2)-1)';
    end

%% Takes the projection relative to the snip keys
    
        % Taking subset of the region
        % forward{i} = convSquare(startH:startH+size(testPoint,1)-1,startL:startL+size(testPoint,2)-1);
        
        % Regular method of convolution.
        % forward{i} = conv2(testPoint,rot90(weights(:,:,1,i),2),'same')+bias(i);
% Reshaping Input
                                                                            % some reason snipKeys need to be transposed cause dim keys are row
                                                                            % vector.
for i = 1: length(weights)
    X1conv1{i,2}  = X1conv1{i,1}(snipKeys{i}',:,:);
end
    % for i = 1:length(net1.Layers(2).Weights)
%     x = [x ;reshape(forward{i},[1 size(forward{i},1)*size(forward{i},2)])];
% end

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

for j = 1: strideVert : size(testPoint,1)
    for i = 1: strideHorz :size(testPoint,2)
        zeroMatrix = zeros(size(testPoint));
        zeroMatrix(j:j+poolHeight-1,i:i+poolLength-1) = zeroMatrix(j:j+poolHeight-1,i:i+poolLength-1)+poolMatrix;

        poolingShape(k,:) = reshape(zeroMatrix',[size(testPoint,1)*size(testPoint,2) 1])';
        k = k + 1;
    end
end

%% Pooling Calculation
for i = 1: length(weights)
        % Create unique keys per conv layer
        %%change pooling key definition hard coded
        poolKeys{i} = genMatIndexKeys(sprintf('P_%d',i),1,1, 14,14);
end
for i = 1: length(weights)
    x{i,1} = X1conv1{i,2}.transform(zeros(196,1),poolingShape,snipKeys{i},poolKeys{i});
end
% sparsePool = kron(sparse(eye(200)),sparse(poolingShape));
% x = x.transform(sparse(zeros(39200,1)),sparsePool,snipKeys,memZono.genKeys('P',1:39200))
% x = poolingShape*x';    % calcs fine
%%
X = xConvOut;
% Taking ReLU with ReLUNN 
% First Relu Layer (57600*50)
% W * x + B
Ws = {double(net1.Layers(4).Weights), double(net1.Layers(6).Weights), double(net1.Layers(8).Weights), double(net1.Layers(10).Weights)};
bs = {double(net1.Layers(4).Bias), double(net1.Layers(6).Bias), double(net1.Layers(8).Bias), double(net1.Layers(10).Bias)};

[NN, Z] = reluNN(X, Ws, bs, 1000);

range = Z.Z('y_1').bounds;
% %% Relu layer 4-5
% weights = net1.Layers(4,1).Weights;
% x = reshape(x,[size(weights,2) 1]);
% bias = net1.Layers(4,1).Bias;
% 
% % linear mapping
% forwardRelu1 = weights*x+bias;
% 
% % relu activation function
% activeRelu1 = max(0,forwardRelu1);
% 
% %% Relu layer 6-7
% weights = net1.Layers(6,1).Weights;
% bias = net1.Layers(6,1).Bias;
% 
% % linear mapping
% forwardRelu2 = weights*activeRelu1 + bias;
% 
% % relu activation function
% activeRelu2 = max(0,forwardRelu2);
% 
% %% Relu layer 8-9
% weights = net1.Layers(8,1).Weights;
% bias = net1.Layers(8,1).Bias;
% 
% % linear mapping
% forwardRelu3 = weights*activeRelu2 + bias;
% 
% % relu activation function
% activeRelu3 = max(0,forwardRelu3);
% 
% %% Relu layer 10
% weights = net1.Layers(10,1).Weights;
% bias = net1.Layers(10,1).Bias;
% 
% % linear mapping
% lastLayer(start) = weights*activeRelu3 + bias;
% % end

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
function keys = genMatIndexKeys(key,Is, Js, I,J)
    keys = {};
    for j = Js:J
        for i = Is: I
        keys = [keys,memZono.genKeys(sprintf("%s_%d",key,i),j)];
        end
    end
end