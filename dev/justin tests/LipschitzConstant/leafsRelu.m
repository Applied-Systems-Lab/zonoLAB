clear; clc; close all;

load('sincos_20_10_10.mat',"net","NN")
load("vertciesAndFacets.mat")

for i = 2:2:length(net.Layers)
    Weights{i/2} = net.Layers(i,1).Weights;
    Bias{i/2} = net.Layers(i,1).Bias;
end

for leaf = 1: length(LCfacets)
    
    allVerticesForLeaf = [];
    vertexNum = size(rmmissing(LCfacets(leaf,:)),2);
    for i = 1 : vertexNum
        allVerticesForLeaf = [allVerticesForLeaf; LCvertices(LCfacets(leaf,i),:)]; 
    end
    
    averageVert(leaf,:) = mean(allVerticesForLeaf);

    layer{leaf,1} = [averageVert(leaf,1:2)'];

    for i = 1:length(Weights)
        layer{leaf,i+1} = Weights{i}*layer{leaf,i} + Bias{i};

        if(i == length(Weights))
            break;
        end
    
        for j = 1: length(layer{leaf,i+1})
            if(layer{leaf,i+1}(j) < 0)
                layer{leaf,i+1}(j) = 0;
            end
        end
    end
end

%% Graphing Layers
f = figure('name',"Plotting Leaves")
tLen = 20;

% [sorted, I] = sort(averageVert);
% Zi = Zi(I)';
% layer = layer(I,:);

%% Graphing Leaves
sp1 = subplot(2,1,1);

optPlot = plotOptions('Display','on','SolverOpts',{},'FaceColor',[1 0 0],'FaceAlpha',0.7);
for leaf = 1 : length(Zi)
    subplot(2,1,1);
    plot(Zi{leaf},optPlot);
    axis([-5 5 -5 5])
    subplot(2,1,2);
    for i = 1 : size(layer,2)
        for j = 1 : length(layer{leaf,i})
            % Add something for layer 1
            if (i == 1)
                p3 = plot(i,((j-1)*tLen/length(layer{leaf,i}))+(tLen/(2*length(layer{leaf,i}))),"o", "MarkerSize",10,'MarkerEdgeColor','k','MarkerFaceColor','k');
            %Add something for last layer
            elseif(i == size(layer,2))
                p4 = plot(i,((j-1)*tLen/length(layer{leaf,i}))+(tLen/(2*length(layer{leaf,i}))),"o", "MarkerSize",10,'MarkerEdgeColor','k','MarkerFaceColor','b');
            else
                % Any Layer
                if(layer{leaf,i}(j) > 0)
                    p1 = plot(i,((j-1)*tLen/length(layer{leaf,i}))+(tLen/(2*length(layer{leaf,i}))),"o", "MarkerSize",10,'MarkerEdgeColor','k','MarkerFaceColor','g');
                else
                    p2 = plot(i,((j-1)*tLen/length(layer{leaf,i}))+(tLen/(2*length(layer{leaf,i}))),"o", "MarkerSize",10,'MarkerEdgeColor','k','MarkerFaceColor','r');
                end
            end

            hold on
        end
    end
    legend([p1, p2], ["On","Off"]);
    axis off

    subplot(2,1,1);
    view(0,90)

    pause(1)
    saveas(gcf, append(sprintf("%d", leaf),".png"))
%     while(~waitforbuttonpress)
%     end
   
%     cla(sp1);
end