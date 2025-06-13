% Programmer Name: Justin Chen
% The following code is a updated version of the 2022-11-29 version of the
% mnist command file form the L4DC papers.

% List of Changes


clear; close all;
clc;

% True trains the mnist data
% False load a NN trained on 784-5-5-1
fromScratch = true;

%% Training Data
load('classification_data.mat')
N = size(data,2);
x1 = class1(1,:);
y1 = class1(2,:);
classification1 = class1(3,:);
classification1(1,:) = 1;

x2 = class2(1,:);
y2 = class2(2,:);
classification2 = class2(3,:);

x1_train = [class1(1,:) class2(1,:)];
x2_train = [class1(2,:) class2(2,:)];

% Training Parameters
input_train = [x1_train;x2_train];
% output_train = [zeros(1,size(classification1,2)) classification2];
output_train = [classification1 zeros(1,size(classification2,2));zeros(1,size(classification1,2)) classification2];

%% Training Network
if(fromScratch)
    layers = [featureInputLayer(2) 
             fullyConnectedLayer(10)
             reluLayer
             fullyConnectedLayer(10)
             reluLayer
             fullyConnectedLayer(10)
             reluLayer
             fullyConnectedLayer(5)
             reluLayer
             fullyConnectedLayer(2)             % Was changed from 1 to 2?
             regressionLayer];
         
    options = trainingOptions("sgdm",'Momentum', 0.95, MaxEpochs = 10000,Plots='training-progress');
         
    net = trainNetwork(input_train', output_train', layers, options);

    save('200dpClass.mat')
else
    % load('2000dp10000epMIMO1+1Classification.mat')
end


%% Generate data set

input_test = data(1:2,1:0.2*N);

output_test = predict(net,input_test');
output_validate = data(3,1:0.2*N)';

%%
% output_validate(find(output_validate(:,1)==-1),1)=0;

diff = output_test - output_validate;
test = abs(output_test);

for i = 1: length(test)
    if(test(i,1)>test(i,2))
        test(i,:)=[1,0];
    elseif (test(i,1)<test(i,2))
        test(i,:)=[0,1];
    else
        test(i,:)=[0,0];
    end
end

accuracy = (1-sum(abs(output_validate-test))/(0.2*N))*100;

%% Hyb zono
domain = [-1,1,-1,1];
Ws = [];
bs = [];

for i = 1: floor(length(layers)/2)
    Ws = [Ws {double(net.Layers(2*i).Weights)}];
    bs = [bs {double(net.Layers(2*i).Bias)}];
end

cdomain = num2cell(domain);
[x1_min, x1_max, x2_min, x2_max] = cdomain{:};
g11 = (x1_max - x1_min)/2;
g22 = (x2_max - x2_min)/2;


Gx = diag([g11, g22]);

cx = zeros(2, 1);
X = hybZono(Gx, [], cx, [], [], []);
X = memZono(X,'X');
a = 1000;

[NN,Y] = reluNN(X,Ws,bs,a);

%%
NNc1 = NN(NN.dimKeys(1:3));
NNc2 = NN([NN.dimKeys(1:2) NN.dimKeys(4)]);

%%
AZ =  -124;
EL =   50;

fig = figure;

f = @(x) 15.*(x-.5).*(x+1).*x.*(x+.5).*sin(1.2*(x-.25))-.35;
fplot(f, [-1 1],'Color',[0.5 0.5 0.5],'LineWidth',3)
hold on
% plot(goods(1,:), goods(2,:), 'r.', 'MarkerSize', 10)
% plot(bads(1,:), bads(2,:), 'b.','MarkerSize',10)

axis([-1 1 -1 1])
box off


for i = 1:N
    if(data(3,i) == 1)
        scatter3(data(1,i),data(2,i),1, 40,'b.');
    elseif(data(3,i) == 0)
        scatter3(data(1,i),data(2,i),0, 40,'r.');
    end
end

%%
plot(NNc2.Z(NNc2.dimKeys),'k',0.2)
view([AZ,EL])

%% Retraining

