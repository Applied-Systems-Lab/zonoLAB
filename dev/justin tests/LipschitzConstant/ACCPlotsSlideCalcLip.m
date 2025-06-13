%% Code of Plots for ACC
clc
clear
close all

AZ =  -37.5000;
EL =   30;

f = figure;
f.Position = [603 346.6000 826.8000 577.4000];

load('6_2linear_4_5_4.mat','NN')
load('Lip6-2Lin.mat','L')

tempZono = NN.Z(NN.dimKeys);
leaves = round(tempZono.getLeaves({}));
[ngb, num_leaves] = size(leaves);
[LCvertices, LCfacets] = plotHybZono3D(tempZono,{})

for leaf = 1 : num_leaves
    Zi{leaf,1} = conZono(tempZono.Gc,tempZono.c+tempZono.Gb*leaves(:,leaf),tempZono.Ac,tempZono.b-tempZono.Ab*leaves(:,leaf));
end

plot(Zi{1},'k',0.2)
hold on
scatter3(LCvertices(1:5,1),LCvertices(1:5,2),LCvertices(1:5,3),45,'or')
title(sprintf('Leaf: %d Slope: %.2f',1,L(1,1)))

saveas(gcf,'findPoints.png')