clear
clc
close all

load('reLUMnistwConv1111111.mat')

str = {'conv', 'relu'};
j = 1;

for i = 1: length(net1.Layers)
    if(contains(net1.Layers(i).Name,str))
        L{j} = net1.Layers(i).Name;
        j = j + 1;
    end
end

disp(L);

%%
j = 1;
for i = 1: length(net1.Layers)
    if contains(net1.Layers(i).Name,str{1})
            Ws{j} = net1.Layers(i).Weights;
            bs{j} = net1.Layers(i).Bias;
            j = j + 1;
    elseif contains(net1.Layers(i).Name,str{2})
            Ws{j} = net1.Layers(i-1).Weights;
            bs{j} = net1.Layers(i-1).Bias;
            j = j + 1;
    end
end

% Gets the last fc layer used for the regression layer.
Ws{j} = net1.Layers(i).Weights;
bs{j} = net1.Layers(i).Bias;

%%
ones_data = cat(4,XValidation(:,:,:,XValidLabels == 1) , XTest(:,:,:,XTestLabels == 1) , XTrain(:,:,:,XTrainLabels == 1));
sevens_data = cat(4,XValidation(:,:,:,XValidLabels == 7) , XTest(:,:,:,XTestLabels == 7) , XTrain(:,:,:,XTrainLabels == 7));
data = cat(4,ones_data,sevens_data);

%%

c = mean(data,4);
M = max(data,[],4);
m = min(data,[],4);
s = std(data,0,4);

c = reshape(c, [28*28,1]);
m = reshape(m, [28*28,1]);
M = reshape(M, [28*28,1]);
Gc = diag((M-m)/2);

Z = hybZono(Gc,[],c,[],[],[]);
Z = memZono(Z,'Z');

%%

a = 1000;       % Large Number

convStep(X,Ws,bs,a,L);
