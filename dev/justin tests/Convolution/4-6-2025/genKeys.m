clear
clc
close all

load('toplitest.mat')

%% Contents passed into the function
% Weights and Bias might be seperate functions
weights = net1.Layers(2,1).Weights;
bias = net1.Layers(2,1).Bias;
input = avg_ones_data;

convHeight = size(weights(:,:,1,1),1);
convLength = size(weights(:,:,1,1),2);

% Creates empty matrix accounting for padding?
Empty = zeros(size(input,1)+2*(convHeight-1),size(input,2)+2*(convLength-1));

%% Condensing for loops on gen keys
inputHeight = 28;
inputWidth = 28;
outKeys = [];
tic
for i = 1: length(weights)
    % Create unique keys per conv layer
    outKeys{i} = genMatIndexKeysRColumnRow(sprintf('C_%d',i),1,1,inputHeight+convHeight-1,inputWidth+convLength-1)';
end
toc
disp('end')

%%
a = length(weights);
b = inputHeight+convHeight-1;
c = inputWidth+convLength-1;
tic

convKey = @(w,i,j) sprintf('%s_%d_%d',w,i,j);
convKeys = @(w,row,col) arrayfun(@(i,j) convKey(w,i,j), row, col , UniformOutput=false);

% convKey = @(w,i,j) sprintf('C_%d_%d_%d',w,i,j);
% convKeys = @(w,I1s,I2s,I1e,I2e) arrayfun(@(i,j) convKey(w,i,j), kron(I1s:I1e,ones(1,I2e-I2s+1)),repmat(I2s:I2e,[1,I1e-I1s+1]) , UniformOutput=false);
for i = 1:a
    outKeys2{i} = convKeys(sprintf('C_%d',i),kron(1:b,ones(1,c-1+1)),repmat(1:c,[1,b-1+1]))'; 
end
out = vertcat(outKeys2{:});
toc

%%
tic
inputKeys = genMatIndexKeysRRowColumn('Z',1,1,inputHeight+2*(convHeight-1),inputWidth+2*(convLength-1))';
toc

tic
inputHwpad = inputHeight+2*(convHeight-1);
inputWwpad = inputWidth+2*(convLength-1);
hStart = 1;
wStart = 1;

inputKeys2 = convKeys('Z',repmat(1:inputWwpad,[1,inputHwpad-hStart+1]),kron(1:inputHwpad,ones(1,inputWwpad-wStart+1)))';
toc

isequal(inputKeys,inputKeys2)

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