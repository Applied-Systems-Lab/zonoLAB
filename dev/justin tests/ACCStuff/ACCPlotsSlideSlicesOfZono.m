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

t = tiledlayout(3,4);
nexttile(t,[3,3])
for leaf = 1: num_leaves
        plot(Zi{leaf},'k',0.2)
end
plot(Zi{3},'g',0.8)
plot(Zi{5},'b',0.8)
plot(Zi{7},'r',0.8)

set(gca,XTick = [-3, 0, 3],YTick = [-3, 0, 3],ZTick = [-20, 0, 20], FontSize = 18);
xlabel('$x_1$','interpreter','latex','FontSize',28)
ylabel('$x_2$','interpreter','latex','FontSize',28)
zlabel('$\mathcal{F}(x)$','interpreter','latex','FontSize',28)

nexttile(t,[4])
plot(Zi{3},'g',0.8)
title(sprintf('Slope: %.2f',L(3,1)))

nexttile(t,[8])
plot(Zi{7},'r',0.8)
title(sprintf('Slope: %.2f',L(7,1)))

nexttile(t,[12])
plot(Zi{5},'b',0.8)
title(sprintf('Slope: %.2f',L(5,1)))

%%
annotation('arrow',[0.5 0.8],[0.4 0.8],'Color', 'k','LineWidth',2);     % first subplot
annotation('arrow',[0.5 0.8],[0.335 0.5],'Color', 'k','LineWidth',2);     % second subplot
annotation('arrow',[0.4 0.8],[0.25 0.25],'Color', 'k','LineWidth',2);     % third subplot

%%
saveas(gcf,'slicesofZono.png')