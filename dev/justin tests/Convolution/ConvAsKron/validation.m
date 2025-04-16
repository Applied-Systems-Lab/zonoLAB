clear
clc
close all

load("poolShape.mat")

% How poolingShape is generated
% %% Average Pooling Layer
% poolLength = net1.Layers(3).PoolSize(1);
% poolHeight = net1.Layers(3).PoolSize(2);
% strideVert = net1.Layers(3).Stride(1);
% strideHorz = net1.Layers(3).Stride(2);
% 
% 
%% Number of cells averaged
pool = 1/(poolLength*poolHeight);

poolMatrix = pool*ones(2,2);

%%
k = 1;
tic
for j = 1: strideVert : size(forward{1},1)
    for i = 1: strideHorz :size(forward{1},2)
        zeroMatrix = zeros(size(forward{i}));
        zeroMatrix(j:j+poolHeight-1,i:i+poolLength-1) = zeroMatrix(j:j+poolHeight-1,i:i+poolLength-1)+poolMatrix;

        poolingShape(k,:) = reshape(zeroMatrix',[size(forward{i},1)*size(forward{i},2) 1])';
        k = k + 1;
    end
end
toc

%% Very basic example. Not implemented dynamically
tic
% The content of the pooling layer specifically the first row of the
% pooling 
A = [0.25 0.25];

% Current guess is this is the output size of the matrix as a square
B = eye(14);

% The number of rows in the pooling layer
C = [1,1];

E = kron(B,kron(C, kron(B,A)));
toc
%%
% x = poolingShape*x';

% x is the stacked vectorized input layer

%% Goal to compare poolingShape with a kron version

% Given A, B: if A is the contents and B is the structure,
% kron(A,B) will give the pooling shape as a single matrix if the size is
% defined correspondingly

A = zeros(1,28);
A(1,1) = 0.25;
A(1,2) = 0.25;

B1 = ones(1,2);

C1 = kron(B1,A);

B2 = eye(14);

C2 = kron(C1,B2);

B3 = eye(14);

C3 = kron(C2,B3);
