% Sorting all the arrays
clear; clc; close all;

load('sincos_20_10_10.mat',"net","NN")
load("vertciesAndFacets.mat")
load("cursor_info.mat")

for i = 2:2:length(net.Layers)
    Weights{i/2} = net.Layers(i,1).Weights;
    Bias{i/2} = net.Layers(i,1).Bias;
end

for leaf = 1: length(cursor_info)
    
    layer{leaf,1} = [cursor_info(leaf).Position'];

    for i = 1:length(Weights)
        layer{leaf,i+1} = Weights{i}*layer{leaf,i} + Bias{i};

        if(i == length(Weights))
            break;
        end
        tempLayer(leaf,i+1) = layer(leaf,i+1);
        for j = 1: length(layer{leaf,i+1})
            if(layer{leaf,i+1}(j) < 0)
                layer{leaf,i+1}(j) = 0;
            end
        end
    end
end

newLay = zeros(leaf,43);

for i = 1: leaf
    newLay(i,:) = [layer{i,1};layer{i,2};layer{i,3};layer{i,4};layer{i,5}]';
end

[newLay, I] = sortrows(newLay,[3:42]);

layer = layer(I,:);
Zi = Zi';
Zi = Zi(I);

% save("leafRELU.mat")

%% Graphing Layers
f = figure('name',"Plotting Leaves","Position",[454 276 1705 838]);

leafAnno =    annotation( 'textbox' , [0.5 0.83 0.1 0.1], ...
                'String' , " ", ...
                'FitBoxToText' , 'on', ...
                'LineStyle' , 'none', ...
                'FontSize' , 16, ...
                'FontWeight' , 'bold');

firstNeuro = annotation('textarrow',[0.624 0.634],[0.12 0.13],'Color', 'k',"String"," First Neuron ");
 lastNeuro = annotation('textarrow',[0.835 0.825],[0.893 0.883],'Color', 'k',"String"," Last Neuron ");
tLen = 20;
% [sorted, I] = sort(averageVert);
% Zi = Zi(I)';
% layer = layer(I,:);

%% Graphing Leaves
t = tiledlayout(1,2);
a = 1; % temp variable for text boxes

optPlot = plotOptions('Display','on','SolverOpts',{},'FaceColor',[1 0 0],'FaceAlpha',0.7);
for leaf = 1 : length(cursor_info)
    nexttile(1)
    plot(Zi{leaf},optPlot);
    hold on
    view(0,90)
    axis([-5 5 -5 5])
    reLuPlot = nexttile(2);
    cla(reLuPlot)
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
    l = legend([p1, p2, p3, p4], ["On","Off", "Inputs", "Outputs"]);
    l.Location = 'southeast';
    
    axis off

    leafAnno.String = sprintf('Leaf: %d',I(leaf));

%     pause(0.01);
%     waitforbuttonpress
%     saveas(gcf,(sprintf('%d_%d.png',leaf,I(leaf))));
end

nexttile(1)
hold on
for i = 1: length(cursor_info)
    plot(newLay(i,1),newLay(i,2),'k.', "MarkerSize",10)
    hold on
end
plot(newLay(67,1),newLay(67,2),'gx', "MarkerSize",10)