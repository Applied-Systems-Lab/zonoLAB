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

for leaf = 1 : num_leaves
    Zi{leaf,1} = conZono(tempZono.Gc,tempZono.c+tempZono.Gb*leaves(:,leaf),tempZono.Ac,tempZono.b-tempZono.Ab*leaves(:,leaf));
end

t = tiledlayout(3,3);
for leaf = 1: num_leaves
    nexttile(t,[leaf])
    plot(Zi{leaf},'k',0.2)
    title(sprintf('Leaf: %d Slope: %.2f',leaf,L(leaf,1)))
end

saveas(gcf,'allLeaves.png')