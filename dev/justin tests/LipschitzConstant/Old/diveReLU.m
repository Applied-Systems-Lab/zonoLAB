clear; clc; close all;

load('sincos_20_10_10.mat',"net","NN")
load("leaf189singleton","vertices")
load("leaf189CursorPoint.mat","cursor_info")
load("vertciesAndFacets.mat")

x = vertices;

prediction = net.predict([x(1),x(2)]);

x(3) - prediction

for i = 2:2:length(net.Layers)
    Weights{i/2} = net.Layers(i,1).Weights;
    Bias{i/2} = net.Layers(i,1).Bias;
end

layer{1} = [x(1);x(2)];

for i = 1:length(Weights)
    layer{i+1} = Weights{i}*layer{i} + Bias{i};

    if(i == length(Weights))
        break;
    end

    for j = 1: length(layer{i+1})
        if(layer{i+1}(j) < 0)
            layer{i+1}(j) = 0;
        end
    end
end

leaf = NN.Z.getLeaves({});
stackedLayer = [layer{2};layer{3};layer{4}];
reLu189 = [stackedLayer, leaf(:,189)];

% %% Calculating the Relu for another point
% 
% for k = 1 : length(cursor_info)
%     x = cursor_info(k).Position;
% 
%     for i = 2:2:length(net.Layers)
%         Weights{i/2} = net.Layers(i,1).Weights;
%         Bias{i/2} = net.Layers(i,1).Bias;
%     end
%     
%     layer{k+1,1} = [x(1);x(2)];
% 
%     for i = 1:length(Weights)
%         layer{k+1,i+1} = Weights{i}*layer{k+1,i} + Bias{i};
% 
%         if(i == length(Weights))
%             break;
%         end
% 
%         for j = 1: length(layer{k+1,i+1})
%             if(layer{k+1,i+1}(j) < 0)
%                 layer{k+1,i+1}(j) = 0;
%             end
%         end
%     end
% end
% allReLU = [];
% for i = 2: length(layer)-1
%     allReLU = [allReLU ; layer{:,i}];
% end

%% Graphing Layers
figure('name',"Plotting Leaves")
tLen = 20;

%% Graphing Leaves
sp1 = subplot(2,1,1);

optPlot = plotOptions('Display','on','SolverOpts',{},'FaceColor',[1 0 0],'FaceAlpha',0.7);
for leaf = 1 : length(Zi)
    subplot(2,1,1);
    plot(Zi{leaf},optPlot);
    axis([-5 5 -5 5 -5 5])
    subplot(2,1,2);
    for i = 1 : length(layer)
        for j = 1 : length(layer{i})
            if(layer{i}(j) > 0)
                p1 = plot(i,((j-1)*tLen/length(layer{i}))+(tLen/(2*length(layer{i}))),"o", "MarkerSize",10,'MarkerEdgeColor','k','MarkerFaceColor','g');
            else
                p2 = plot(i,((j-1)*tLen/length(layer{i}))+(tLen/(2*length(layer{i}))),"o", "MarkerSize",10,'MarkerEdgeColor','k','MarkerFaceColor','r');
            end
            hold on
        end
    end
    legend([p1, p2], ["On","Off"]);
    axis off

    subplot(2,1,1);
    view(0,90)
    while(~waitforbuttonpress)
    end
   
%     cla(sp1);
end


