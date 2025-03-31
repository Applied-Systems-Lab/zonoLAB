% Sorting all the arrays
clear; clc; close all;

load("leafRELU.mat")
newLay = zeros(795,43)
for i = 1: leaf
    newLay(i,:) = [layer{i,1};layer{i,2};layer{i,3};layer{i,4};layer{i,5}]';
end

[sorted, I] = sortrows(newLay,[3:42]);

layer = layer(I,:);
Zi = Zi';
Zi = Zi(I);

%% Graphing Layers
f = figure('name',"Plotting Leaves","Position",[454 276 1705 838])

leafAnno =    annotation( 'textbox' , [0.5 0.83 0.1 0.1], ...
                'String' , " ", ...
                'FitBoxToText' , 'on', ...
                'LineStyle' , 'none', ...
                'FontSize' , 16, ...
                'FontWeight' , 'bold');
tLen = 20;
% [sorted, I] = sort(averageVert);
% Zi = Zi(I)';
% layer = layer(I,:);

%% Graphing Leaves
t = tiledlayout(1,2);
a = 1; % temp variable for text boxes

optPlot = plotOptions('Display','on','SolverOpts',{},'FaceColor',[1 0 0],'FaceAlpha',0.7);
for leaf = 1 : 100%length(Zi)
    nexttile(1)
    plot(Zi{leaf},optPlot);
    hold on
    view(0,90)
    axis([-5 5 -5 5])
    nexttile(2)
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
    legend([p1, p2, p3, p4], ["On","Off", "Inputs", "Outputs"]);
    axis off

%     M{leaf} = getframe().cdata;
%     while(~waitforbuttonpress)
%     end
    leafAnno.String = sprintf('Leaf: %d',I(leaf));
    pause(0.01);
%     cla(sp1);
end
% gifFile = 'myAnimation.gif';
% imwrite(M{1}, gifFile, 'DelayTime', 0.5, 'LoopCount', inf);
% 
% for n = 2:numel(M)
%     imwrite(M{n}, gifFile, 'DelayTime', 0.5, 'LoopCount', inf);
% end
