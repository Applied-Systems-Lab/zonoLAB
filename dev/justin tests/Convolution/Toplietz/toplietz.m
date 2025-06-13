clear
clc
close all

%%
load("toplietz.mat")

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