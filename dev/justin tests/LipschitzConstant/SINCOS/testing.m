
clear all; close all; clc;

AZ =  -14.7764;
EL =   -3.2364;

f = figure;
f.Position = [603 347 1301 577];
% load('Papersincos_20_10_10.mat',"NN")
load('Paperlinear_4_5_4.mat',"NN")
subplot(1,2,1)

plot(NN.Z,'r',0.8);
view(AZ,EL)
xlabel('$X_1$','interpreter','latex')
ylabel('$X_2$','interpreter','latex')
zlabel('$output$','interpreter','latex')

load('PaperDiffEpLinear2_4_5_4.mat',"NN")

Zi{leaf,1} = conZono(tempZono.Gc,tempZono.c+tempZono.Gb*leaves(:,leaf),tempZono.Ac,tempZono.b-tempZono.Ab*leaves(:,leaf));
subplot(1,2,2)

plot(NN.Z,'r',0.8);
view(AZ,EL)
xlabel('$X_1$','interpreter','latex')
ylabel('$X_2$','interpreter','latex')
zlabel('$output$','interpreter','latex')

saveas(gcf,'diffEpoch.png')
