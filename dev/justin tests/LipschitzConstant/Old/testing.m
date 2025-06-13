
clear, close all, clc

load('sincos_20_10_10.mat',"net")



%% Problem Points
% load("leaf549.mat")
% load("leaf548.mat")
load("leaf189singleton","vertices")

% x = LCvertices(2187,1:2);     % Leaf 549
% x1 = LCvertices(2183,1:2);     % Leaf 548 1
% x2 = LCvertices(2184,1:2);     % Leaf 548 2
x = vertices(1,1:2);     % Leaf 189

xspan = x(1)-0.0002 : 0.00001: x(1)+0.0002;
yspan = x(2)-0.0002 : 0.00001: x(2)+0.0002;

[X,Y] = meshgrid(xspan,yspan);
Z = arrayfun(@(x,y)net.predict([x,y]),X,Y);

% y = zeros(length(xspan)*length(yspan), 3);
% for i = 1: length(xspan)
%     for j = 1: length(yspan)
%         y(((i-1)*length(yspan))+j,:) = [xspan(i), yspan(j), net.predict([xspan(i), yspan(j)])];
%     end
% end

figure()
surf(X,Y,Z);
hold on
plot3(x(1,1),x(1,2),net.predict([x(1,1),x(1,2)]),'ro');
hold on
plot3(vertices(1,1),vertices(1,2),vertices(1,3),'kx');
hold off

%%
load('Papersincos_20_10_10.mat',"NN")
[sincosSlope, mNorm]= calcLC(NN,1e-7);
avgSlope = mean(sincosSlope)
stdSinCos = std(sincosSlope(:,1),1)

load('Papersincos2X_20_10_10.mat',"NN")
[sincos2XSlope, m2Norm] = calcLC(NN,1e-7);
avgSlope2X = mean(sincos2XSlope)
stdSinCos2X = std(sincos2XSlope(:,1),1)

%%

% input_test are the same for both
load('Papersincos_20_10_10.mat')
sincosInput_Test = input_test;
sincosOutput_test = predict(net,input_test);

load('Papersincos2X_20_10_10.mat')
sincos2XInput_Test = input_test;
sincos2XOutput_test = predict(net,input_test);

% theoretically output test should be 2x1 if we take the difference between
% the two functions.

% cos(x1) + sin(x2) = h1
% cos(x1) + sin(x2) + 2x2 = h2
% h2 - h1 = 2x2
% or double the x2 domain input.

diff = sincos2XOutput_test - 2*input_test(:,2); % should be 2x2
error = sincosOutput_test - diff;