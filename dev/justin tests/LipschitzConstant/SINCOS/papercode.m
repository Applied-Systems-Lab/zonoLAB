%% Section for producting plots for sin and cos example
clc
clear
close all

AZ =  -37.5000;
EL =   30;

f = figure;
f.Position = [603 347 1301 577];

load('Papersincos_20_10_10.mat')

subplot(1,2,1)
plot(NN.Z,'r',0.8)
set(gca,XTick = [-3, 0, 3],YTick = [-3, 0, 3],ZTick = [-2, 0, 2], FontSize = 18);
xlabel('$x_1$','interpreter','latex','FontSize',28)
ylabel('$x_2$','interpreter','latex','FontSize',28)
zlabel('$\mathcal{F}(x)$','interpreter','latex','FontSize',28)

load('Papersincos2X_20_10_10.mat')
subplot(1,2,2)
plot(NN.Z,'r',0.8)
set(gca,XTick = [-3, 0, 3],YTick = [-3, 0, 3],ZTick = [-10, 0, 10], FontSize = 18);
xlabel('$x_1$','interpreter','latex','FontSize',28)
ylabel('$x_2$','interpreter','latex','FontSize',28)
zlabel('$\mathcal{F}(x)$','interpreter','latex','FontSize',28)

% Note there might be statistical significance that there is less leaves in
% the 2x2 version

%% Section for clarifying hybzono and facets
clc
clear
close all

AZ =  -14.7764;
EL =   -3.2364;

f = figure;
f.Position = [603 347 1301 577];
load('Paperlinear_4_5_4.mat')
tempZono = NN.Z;
leaves = round(tempZono.getLeaves({}));
[ngb, num_leaves] = size(leaves);

for leaf = 1 : num_leaves
    Zi{leaf,1} = conZono(tempZono.Gc,tempZono.c+tempZono.Gb*leaves(:,leaf),tempZono.Ac,tempZono.b-tempZono.Ab*leaves(:,leaf));
end

subplot(1,2,1)
for leaf = 1: num_leaves
    if any(leaf==[3])
        continue
    end
        plot(Zi{leaf},'r',0.8)
end

plot(Zi{3},'r',0.8)
set(gca,XTick = [-3, 0, 3],YTick = [-3, 0, 3],ZTick = [-20, 0, 20], FontSize = 18);
xlabel('$x_1$','interpreter','latex','FontSize',28)
ylabel('$x_2$','interpreter','latex','FontSize',28)
zlabel('$\mathcal{F}(x)$','interpreter','latex','FontSize',28)

load('PaperDiffEpLinear2_4_5_4.mat')
tempZono = NN.Z;
leaves = round(tempZono.getLeaves({}));
[ngb, num_leaves] = size(leaves);

for leaf = 1 : num_leaves
    Zi{leaf,1} = conZono(tempZono.Gc,tempZono.c+tempZono.Gb*leaves(:,leaf),tempZono.Ac,tempZono.b-tempZono.Ab*leaves(:,leaf));
end

subplot(1,2,2)
for leaf = 1: num_leaves
    if any(leaf==[16])
        continue
    end
        plot(Zi{leaf},'r',0.8)
end
plot(Zi{16},'r',0.8)
set(gca,XTick = [-3, 0, 3],YTick = [-3, 0, 3],ZTick = [-20, 0, 20], FontSize = 18);
xlabel('$x_1$','interpreter','latex','FontSize',28)
ylabel('$x_2$','interpreter','latex','FontSize',28)
zlabel('$\mathcal{F}(x)$','interpreter','latex','FontSize',28)

%% Linear case with zoom in leaf

clc
clear
close all

AZ = -37.5000;

EL = 30;

f = figure;
f.Position = [0 0 1000 1000];

% color = jet(3);
color = eye(3);

% load('Paperlinear_4_5_4.mat')
load('PaperDiffEpLinear10_4_5_4.mat')

tempZono = NN.Z;
leaves = round(tempZono.getLeaves({}));
[ngb, num_leaves] = size(leaves);

for leaf = 1 : num_leaves
    Zi{leaf,1} = conZono(tempZono.Gc,tempZono.c+tempZono.Gb*leaves(:,leaf),tempZono.Ac,tempZono.b-tempZono.Ab*leaves(:,leaf));
end

t = tiledlayout(4,3);
nexttile(t,[3,3])
hold on

r = 1;
g = 3;
b = 9;

for leaf = 1: num_leaves
    if any(leaf==[r,g,b])
        continue
    end
        plot(Zi{leaf},'k',0.1)
end

plot(Zi{r},color(1,:),0.7)
plot(Zi{g},color(2,:),0.7)
plot(Zi{b},color(3,:),0.7)
set(gca,XTick = [-4, 0, 4],YTick = [-4, 0, 4],ZTick = [-20, 0, 20], FontSize = 18);
xlabel('$x_1$','interpreter','latex','FontSize',28)
ylabel('$x_2$','interpreter','latex','FontSize',28)
zlabel('$F(x)$','interpreter','latex','FontSize',28)

nexttile
plot(Zi{g},color(2,:),0.7)
set(gca,XTick = [],YTick = [],ZTick = []);

nexttile
plot(Zi{r},color(1,:),0.7)
set(gca,XTick = [],YTick = [],ZTick = []);

nexttile
plot(Zi{b},color(3,:),0.7)
set(gca,XTick = [],YTick = [],ZTick = []);

annotation('arrow',[0.27 0.25],[0.59 0.25],'Color', 'k','LineWidth',3);     % first subplot
annotation('arrow',[0.475 0.47],[0.405 0.25],'Color', 'k','LineWidth',3);     % second subplot
annotation('arrow',[0.65 0.75],[0.65 0.25],'Color', 'k','LineWidth',3);     % third subplot

%% Tranlateing back to sin cos
%% Display f based on data points
% Similar to normal graphing methods
figure('Name',name)
subplot(1,4,1)
surf(X1_train, X2_train, Y_train-2*X2_train, 'EdgeColor', 'none');
title('Actual Function')

%% Validate the Network
output_test = predict(net,input_test); 
Y_test = reshape(output_test,p,q);

% The approximation calculated from the training
subplot(1,4,2)
surf(X1_test,X2_test,Y_test-2*X2_test, 'EdgeColor', 'none');
title('Neural Net Approx')

% Difference between the first and second plot
subplot(1,4,3)
surf(X1_test,X2_test,Y_test-Y_validate, 'EdgeColor', 'none');
title('Approximation Error')

%% Calculate Average Error
% Observe how bad the error is across the entire mapping
avgEr = rmse(Y_test,Y_validate,'all');
fprintf('RMSE: %d\n',avgEr);

%% Plot neural network output space 

[x1_min, x1_max] = deal(double(min(x1_test)), double(max(x1_test)));
[x2_min, x2_max] = deal(double(min(x2_test)), double(max(x2_test)));
domain = [x1_min,x1_max,x2_min,x2_max];
Ws = [];
bs = [];

for i = 1: floor(length(layers)/2)
    Ws = [Ws {double(net.Layers(2*i).Weights)}];
    bs = [bs {double(net.Layers(2*i).Bias)}];
end

%% Construct Hybrid Zonotope
tic

% Construct input zonotope
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

fprintf('Zonotope model: ')
toc

%% Plot Hybrid Zonotope

% temp = zono([1 0; 0 1; 1 2],[0;0;0]);
M = eye(3); M(3,2) = -2;
tempZono = NN.Z;

tempM = M*NN.Z;
subplot(1,4,4)
plot(tempM,'r',0.8);
grid on;
title('Hybrid Zonotope')
toc

%% Relative Slope account for area
clc
clear
close all

load('6_2sincos2x_20_10_10.mat')
% sincosSlope = calcLCTwo(NN,1e-7);

% indx = find(sincosSlope == max(sincosSlope(:,4)));

tempZono = NN.Z(NN.dimKeys);
leaves = round(tempZono.getLeaves({}));
[ngb, num_leaves] = size(leaves);
for leaf = 1 : num_leaves
    Zi{leaf,1} = conZono(tempZono.Gc,tempZono.c+tempZono.Gb*leaves(:,leaf),tempZono.Ac,tempZono.b-tempZono.Ab*leaves(:,leaf));
end

figure
plot(NN.Z,'k',0.1)

figure
plot(Zi{mod(indx,num_leaves)},'r',0.8)

%%
% 2x version
load('Papersincos2X_20_10_10.mat')
sincosSlope = calcLCTwo(NN,1e-7);

indx = find(sincosSlope == max(sincosSlope(:,4)));

tempZono = NN.Z;
leaves = round(tempZono.getLeaves({}));
[ngb, num_leaves] = size(leaves);
for leaf = 1 : num_leaves
    Zi{leaf,1} = conZono(tempZono.Gc,tempZono.c+tempZono.Gb*leaves(:,leaf),tempZono.Ac,tempZono.b-tempZono.Ab*leaves(:,leaf));
end

figure
plot(NN.Z,'k',0.1)

figure
plot(Zi{mod(indx,num_leaves)},'r',0.8)

%%
%% Plot neural network output space 

load('Paperlinear_4_5_4.mat')

[x1_min, x1_max] = deal(double(min(x1_test)), double(max(x1_test)));
[x2_min, x2_max] = deal(double(min(x2_test)), double(max(x2_test)));
domain = [x1_min,x1_max,x2_min,x2_max];
Ws = [{[1,1;1,1]},{[1;1]},{[1]}];
bs = [{[1,1]},{[1;1]},{[1]}];

tic

% Construct input zonotope
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

fprintf('Zonotope model: ')
toc

tempZono = NN.Z;

figure
plot(tempZono,'r',0.8);
grid on;
title('Hybrid Zonotope')
toc

%% Paper code for classification with Jonah

figure;

f = gcf;
f.Position = [603 347 2068 612];

file = [{'HZ_flat.mat'} {'2000dp10000epClassification.mat'}];
for iterate = 1:length(file)
    subplot(1,2,iterate)
    if iterate == 2
        load('set_up_data_2000.mat')
    else
        load('set_up_data_200.mat')
    end

    hold on
    plot(class1(1,:), class1(2,:), 'r.', 'MarkerSize', 10)
    plot(class2(1,:), class2(2,:), 'b.','MarkerSize',10)
    fp = fplot(f, [-1 1]);
    fp.Color = [0.5 0.5 0.5];
    fp.LineWidth = 3;
    axis([-1 1 -1 1])

    load(file{iterate});

    if(iterate == 1)
        NN = memZono(HZ_flat,'X');
    end

    plot(NN.Z,'k',0.2);
    
    Gc1 = [  0.15 0;
            0 0.25
          ];
    c1 = [0.85;0];
    
    domainClass2 = conZono(Gc1,c1);
    
    plot(domainClass2,'b',0.2)
    
    intersectedClass2 = NN.Z.and(domainClass2,[1 0 0; 0 1 0]);

    Gc1 = [  0.25 0;
        0 0.15
      ];
    c1 = [0;0.85];
    
    domainClass1 = conZono(Gc1,c1);
    
    plot(domainClass1,'r',0.2)
    
    intersectedClass1 = NN.Z.and(domainClass1,[1 0 0; 0 1 0]);
    
    plottingintersectedClass2 = intersectedClass2;
    plottingintersectedClass2.c = plottingintersectedClass2.c+[0;0;0.003];
    plottingintersectedClass1 = intersectedClass1;
    plottingintersectedClass1.c = plottingintersectedClass1.c+[0;0;0.01];
    
    plot(plottingintersectedClass2,'b',0.5)
    plot(plottingintersectedClass1,'r',0.5)
    
    xlabel('$x_1$','interpreter','latex','FontSize',28)
    ylabel('$x_2$','interpreter','latex','FontSize',28)
    zlabel('$\mathcal{F}(x)$','interpreter','latex','FontSize',28)
        
    AZ =    -123.7833;
    EL =   36.5541;
    
    view(AZ,EL);
end

%%
annotation('textbox',[0.09,.7, .2,.2],'String','(a)','LineStyle','none','interpreter','latex','FontSize',20);
annotation('textbox',[.367,.7, .2,.2],'String','(b)','LineStyle','none','interpreter','latex','FontSize',20);
annotation('textbox',[.65,.7, .2,.2],'String','(c)','LineStyle','none','interpreter','latex','FontSize',20);

%%
clear
clc
close all

load('mimoNN.mat')

NNc1 = NN(NN.dimKeys(1:3));
NNc2 = NN([NN.dimKeys(1:2) NN.dimKeys(4)]);

figure
plot(NNc1.Z,'b', 0.1)
plot(NNc2.Z,'r', 0.1)

%%
Gc1 = [  0.15 0;
        0 0.25
      ];
c1 = [0.85;0];

domainClass2 = conZono(Gc1,c1);

plot(domainClass2,'b',0.2)

intersectedClass2 = NNc2.Z.and(domainClass2,[1 0 0; 0 1 0]);

Gc1 = [  0.25 0;
        0 0.15
      ];
c1 = [0;0.85];

domainClass1 = conZono(Gc1,c1);

plot(domainClass1,'r',0.2)

intersectedClass1 = NNc2.Z.and(domainClass1,[1 0 0; 0 1 0]);

plottingintersectedClass2 = intersectedClass2;
plottingintersectedClass2.c = plottingintersectedClass2.c+[0;0;0.003];
plottingintersectedClass1 = intersectedClass1;
plottingintersectedClass1.c = plottingintersectedClass1.c+[0;0;0.01];

plot(plottingintersectedClass2,'b',0.5)
plot(plottingintersectedClass1,'r',0.5)