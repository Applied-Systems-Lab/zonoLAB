clear
clc
close all

%%
% Theoretically a function:
% Will generate a single matrix that a convolutional function will use

% Currently want to pass in
% 4D array of weights
% Input data size (assumes its consistent across all)



% returns 2D matrix of stacked weights in a toplitez form,
    % matrix convolution as a single matrix operation

load('toplitest.mat')


%% Contents passed into the function
% Weights and Bias might be seperate functions
weights = net1.Layers(2,1).Weights;
bias = net1.Layers(2,1).Bias;
input = avg_ones_data;

%%
% Add in a check to see if all layers are the same
% Assumes all convolutional layers are the same if not need to apply an
% iterate

% Gets the kernel size
convHeight = size(weights(:,:,1,1),1);
convLength = size(weights(:,:,1,1),2);

% Creates empty matrix accounting for padding?
Empty = zeros(size(input,1)+2*(convHeight-1),size(input,2)+2*(convLength-1));

% Corresponds to Lines 68-78
%% WIP: Fixed to be a few lines using topletiz
stackedConvWeights = [];
stackedConvBias = [];
tic
% O(n^3) run time
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
    stackedConvWeights  = [stackedConvWeights; rot90(convShape,2)];
end
toc
disp('created conv shape')

% Compare stackedConv- with stacked-
    % checks out (mean(stackedBias-stackedConvBias))
stackedConvBias = kron(reshape(bias,20,1),ones(1024,1));

% Notes:
% convShape will return a 2D matrix that contains the information of the
% input (columns: 1296 {36 x 36} elements) and the output (rows: 1024 {32 x 32} elements)

%% Final output
    % needs to fix a few hard coded things
% Creates square matrix of the convolutional form
% Still requires snip to the proper output size.
tic
for w = 1: size(weights,4)
    tempw = weights(:,:,:,w);
    % 5 is relative to rows of kernel
    for i = 1:5
        % in actual implementation values should be fixed
            % the zeros padded (31) corresponds to input layer - size of kernel
        temp = triu(toeplitz([tempw(i,:) zeros(1,31)]));
            % 36 is the input layer
            % 32 is output
        temp = kron(diag(ones(36,1),i-1), double(temp(1:32,1:36)));
        out(:,:,i) = temp(1:1024,1:1296);
    end
    convShaped{w} = sum(out,3);
end
convShapeOut = vertcat(convShaped{:});

%%
test = toplizify(weights,Empty);

testaprrox = isapprox(convShapeOut,test);