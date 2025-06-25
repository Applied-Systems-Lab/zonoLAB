clear
clc
close all

%%
load("lipschitzNN.mat")

out1 = lipNN{1,1}.Z('y_1');
stackedNN1 = memoryCartProd(tempX1sigma{1,1},memZono(out1,'y1'));
out7 = lipNN{1,2}.Z('y_1');

stackedNN7 = memoryCartProd(tempX1sigma{1,2},memZono(out7,'y_1'));

%%
load('lipschitzminstV2.mat')
temp = mean(ones_data,2)

find(temp == 0)