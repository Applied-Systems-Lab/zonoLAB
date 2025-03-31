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
start=1;
% for j = 1:2163
    testPoint = reshape(XTest(:,:,1,start),[28 28])-net1.Layers(1).Mean;
    %% For the 7
    weights = net1.Layers(2,1).Weights;
    bias = net1.Layers(2,1).Bias;
    
    % figure
    % hold on
    
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


    % returns a 28 by 28 matrix
    for i = 1: length(weights)
        % subplot(10,10,i)
        
        convLength = size(weights(:,:,1,i),2);
        convHeight = size(weights(:,:,1,i),1);
        Empty = zeros(size(testPoint,1)+2*(convHeight-1),size(testPoint,2)+2*(convLength-1));
        k = 1;
        
        for updown = 1: size(Empty,1)-convHeight+1
            for leftright = 1: size(Empty,2)-convLength+1
                temp = Empty;
                temp(updown:updown+(convHeight-1),leftright:leftright+(convLength-1))= rot90(weights(:,:,1,i),2);
    
                convShape(k,:) = reshape(temp',[size(temp,1)*size(temp,2) 1]);
                k = k + 1;
            end
        end
        
        Ashape = Empty;
        Ashape(convHeight:size(Ashape,1)-convHeight+1,convLength:size(Ashape,2)-convLength+1) = testPoint;
        Ashape = reshape(Ashape',[size(Ashape,1)*size(Ashape,2) 1]);

        % convOut = convShape*Ashape+bias(:,:,i);

        % Full convolution
        convSquare = reshape(rot90(convShape,2)*Ashape+bias(:,:,i),[size(testPoint,2)+convLength-1 size(testPoint,1)+convHeight-1])';

        % Same Convolution partitioning:
        startH = floor(convHeight/2) + 1;
        startL = floor(convLength/2) + 1;
        TESTER = convSquare(startH:startH+size(testPoint,1)-1,startL:startL+size(testPoint,2)-1);

        % forward{i} = conv2(m7,reshape(weights(:,:,1,i),[5,5]),'valid') + bias(i);
        forward{i} = conv2(testPoint,rot90(weights(:,:,1,i),2),'same')+bias(i);
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
    lastLayer(start) = weights*activeRelu3 + bias;
% end

%% Similarity Test

% % FOR ZONOTOPE CASE:
% % mnist.m
% %   lines 141-154           % partitions the data, only need one for now
% %   lines 185 - 230         % finding range
% %   lines 232 - 250 fuzziness
% %   lines 272 - 300
% %   function 320 - 332
% %   function 342 - 355
% %   Shade.m mabye?
% 
% %% Converting NN to zonotopes
% 
% % From the original paper, interval_hull_of_points constructs a zono obj so
% % it does not need to be changed.
% 
% %% Define Input Sets
% X = interval_hull_of_points([XTrain; XTest], 'all');
% 
% ones_data = [ input_train(:,training.labels==1) , input_test(:,test.labels==1) ];
% sevens_data = [ input_train(:,training.labels==7) , input_test(:,test.labels==7) ];
% X_ones = interval_hull_of_points( ones_data, 'all' ); % interval Hull of all ones images
% X_sevens = interval_hull_of_points( sevens_data, 'all' ); % interval Hull of all sevens images
% 
% obvious_ones_data = keep_closest(ones_data,250);
% obvious_sevens_data = keep_closest(sevens_data,250);
% X_obvious_ones = interval_hull_of_points( obvious_ones_data, 1 );
% X_obvious_sevens = interval_hull_of_points( obvious_sevens_data, 1 );
% X_obvious_ones_good = interval_hull_of_points( obvious_ones_data, 0.5 ); %0.1
% X_obvious_sevens_good = interval_hull_of_points( obvious_sevens_data, 0.5 ); %0.1
% 
% %% Average of set data
% figure('Name','Overlaided inputs')
% 
% % Collage of all inputs
% subplot(2,3,1)
% imshow(reshape( X.c , [28,28] ))
% title('all input set')
% 
% % Collage of all 1s
% subplot(2,3,2)
% imshow(reshape( X_ones.c , [28,28] ))
% title('all ones')
% subplot(2,3,3)
% 
% % Collage of the obvious 1s
% imshow(reshape( X_obvious_ones.c , [28,28] ))
% title('obvious ones')
% 
% % Collage of all 7s
% subplot(2,3,4)
% imshow(reshape( X_sevens.c , [28,28] ))
% title('all sevens')
% 
% 
% % Collage of the obvious 7s
% subplot(2,3,5)
% imshow(reshape( X_obvious_sevens.c , [28,28] ))
% title('obvious sevens')
% 
% %% Construct Hybrid Zonotopes
% % Ws and bs are hardcoded and might need to change if the layers are
% % manipulated
% Ws = {double(net.Layers(2).Weights), double(net.Layers(4).Weights), double(net.Layers(6).Weights)};
% bs = {double(net.Layers(2).Bias), double(net.Layers(4).Bias), double(net.Layers(6).Bias)};
% %%
% fprintf('All 1s and 7s:\n')
% Z_range = output_bound(X,Ws,bs);
% fprintf('Range is [ %1.4f , %1.4f ]\n\n', Z_range);
% %%
% fprintf('All 1s:\n')
% Z_ones_range = output_bound(X_ones,Ws,bs);
% fprintf('Range is [ %1.4f , %1.4f ]\n\n', Z_ones_range);
% 
% fprintf('All 7s:\n')
% Z_sevens_range = output_bound(X_sevens,Ws,bs);
% fprintf('Range is [ %1.4f , %1.4f ]\n\n', Z_sevens_range);
% 
% fprintf('All Obvious 1s:\n')
% Z_obvious_ones_range = output_bound(X_obvious_ones,Ws,bs);
% fprintf('Range is [ %1.4f , %1.4f ]\n\n', Z_obvious_ones_range);
% 
% fprintf('All Obvious 7s:\n')
% Z_obvious_sevens_range = output_bound(X_obvious_sevens,Ws,bs);
% fprintf('Range is [ %1.4f , %1.4f ]\n\n', Z_obvious_sevens_range);
% %%
% fprintf('All Obvious (Good) 1s:\n')
% Z_obvious_ones_good_range = output_bound(X_obvious_ones_good,Ws,bs);
% fprintf('Range is [ %1.4f , %1.4f ]\n\n', Z_obvious_ones_good_range);
% 
% fprintf('All Obvious (Good) 7s:\n')
% Z_obvious_sevens_good_range = output_bound(X_obvious_sevens_good,Ws,bs);
% fprintf('Range is [ %1.4f , %1.4f ]\n\n', Z_obvious_sevens_good_range);
% 
% %% 
% subplot(2,3,6)
% plot( [1,1],Z_ones_range, 'r', 'LineWidth', 3); hold on;
% plot( [1,1],Z_obvious_ones_range, 'b', 'LineWidth', 2);
% plot( [7,7],Z_sevens_range, 'm', 'LineWidth', 3);
% plot( [7,7],Z_obvious_sevens_range, 'c', 'LineWidth', 2);
% plot( [4,4],Z_range, 'k', 'LineWidth', 6); hold off;
% legend('1s','Obvious 1s','7s','Obvious 7s','1s and 7s')
% xticks([1,4,7])
% xticklabels({'1s','1s and 7s','7s'})
% xlim([0,8])
% ylabel('Deviation from normal')
% 
% %% Fuzziness from generators
% 
% figure('Name','Adding in uncertainty?')
% subplot(2,3,1)
% imshow(reshape( X.c + X.Gc*ones(size(X.Gc,2),1) , [28,28] ))
% title('all input set')
% subplot(2,3,2)
% imshow(reshape( X_ones.c + X_ones.Gc*ones(size(X_ones.Gc,2),1) , [28,28] ))
% title('all ones')
% subplot(2,3,3)
% imshow(reshape( X_obvious_ones.c + X_obvious_ones.Gc*ones(size(X_obvious_ones.Gc,2),1) , [28,28] ))
% title('obvious ones')
% subplot(2,3,4)
% imshow(reshape( X_sevens.c + X_sevens.Gc*ones(size(X_sevens.Gc,2),1) , [28,28] ))
% title('all sevens')
% subplot(2,3,5)
% imshow(reshape( X_obvious_sevens.c + X_obvious_sevens.Gc*ones(size(X_obvious_sevens.Gc,2),1) , [28,28] ))
% title('obvious sevens')
% 
% %% Final Plots
% 
% figure('Name','Final Plot')
% boxchart(1+0*output_of_1s, output_of_1s, 'BoxEdgeColor','r','MarkerColor','r');
% hold on;
% boxchart(7+0*output_of_7s, output_of_7s, 'BoxEdgeColor','b','MarkerColor','b');
% plot([0.5,7.5],[0,0],'k--');
% hold off;
% 
% xticks([1,7]);
% xticklabels({'1s (test)','7s (test)'})
% yticks([-1,0,1]);
% 
% inset_axes = axes('Position', [0.34 0.14 0.35 0.35]);
% imshow(worst_1)
% pbaspect([1 1 1])
% 
% inset_axes = axes('Position', [0.34 0.55 0.35 0.35]);
% imshow(worst_7)
% pbaspect([1 1 1])
% 
% %%
% X1sigma = {};
% X7sigma = {};
% sigmas = linspace(0.1,0.5,21);
% X1sigma_range = zeros(2,length(sigmas));
% X7sigma_range = zeros(2,length(sigmas));
% for i = 1:length(sigmas)
%     sigma = sigmas(i)
%     X1sigma{i} = interval_hull_of_points( ones_data, sigma );
%     X7sigma{i} = interval_hull_of_points( sevens_data, sigma );
% 
%     X1sigma_range(:,i) = output_bound(X1sigma{i},Ws,bs);
%     X7sigma_range(:,i) = output_bound(X7sigma{i},Ws,bs);
% end
% 
% %%
% figure('Name','Output Range')
% shade(sigmas,X1sigma_range(1,:),'r',sigmas,X1sigma_range(2,:),'r','FillType',[1 2;2 1]);
% hold on;
% shade(sigmas,X7sigma_range(1,:),'b',sigmas,X7sigma_range(2,:),'b','FillType',[1 2;2 1]);
% legend('','','1s','','','7s')
% grid on;
% fontsize(gcf,scale=1.5)
% xlabel('fraction of σ')
% ylabel('prediction interval')
% xlim([min(sigmas),max(sigmas)])
% ylim([-1.4,1.4])
% yticks([-1,0,1])
% 
% %%
% figure()
% subplot(2,1,1)
% imshow(plot_image_versions(X1sigma,5))
% subplot(2,1,2)
% imshow(plot_image_versions(X7sigma,5))
% xlabel('σ')
% ylabel('ξ=1, ξ=-1, ξ∈\{-1,1\}^{784}')
% 
% %%
% figure()
% xi = 2*randi([0,1],size(X.Gc,2),1)-1;
% for i=1:5
% subplot(1,5,i)
% Xsi = X7sigma{abs(sigmas-i/10)<0.0001};
% imshow( reshape( Xsi.c + Xsi.Gc*xi , [28,28] ) )
% end
% 
% % computes the interval hull of a set of points 
% function [Z] = interval_hull_of_points(data,threshold)
%     c = mean(data,2);
%     M = max(data,[],2);
%     m = min(data,[],2);
%     s = std(data,0,2);
%     if threshold == 'all'
%         Gc = diag((M-m)/2);
%     else
%         Gc = diag(threshold*s);
%     end
%     Z = hybZono(Gc,[],c,[],[],[]);
%     Z = memZono(Z,'Z');
% end
% 
% function [range] = output_bound(X,Ws,bs)
%     tic
%     [NN, Z] = reluNN(X, Ws, bs, 1000);
%     toc
%     range = zeros(2,1);
%     % For the specific example, dim order does not matter but could
%     % possibly be an issue
%     bb = Z.Z.boundingBox;
%     ub = bb.c + sum(abs(bb.G),2);
%     lb = bb.c - sum(abs(bb.G),2);
%     range = [lb;ub];
%     % [range]  = Z.Z.bounds;
%     toc
% end