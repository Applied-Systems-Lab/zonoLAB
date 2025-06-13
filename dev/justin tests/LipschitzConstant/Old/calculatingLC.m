% Programmer Name: Justin Chen
% The following program is an attempt to develop a method to calculate the
% lipschtiz constant from trained NN data

clear; close all;
clc;

%% Loading last saved data
fromScratch = true;

%% Loading Data
% Loads pretrained data for 2d and 3d cases

% example = 'hyperp';     % 4d  f = 4*X1 + 2*X2 + 5*X3
example = 'sincos';     % 3d  f = cos(X1) + sin(X2)
% example = 'linear';     % 2d    f = 3*X1

if(example == 'hyperp')
    load('hyperplane_4_4_4.mat',"NN")
elseif(example == 'sincos')
    load('sincos_20_10_10.mat',"NN")
else
    load('linear_4_5_4.mat',"NN")
end

% Current zonoLab version does not have the resolved file, to get around
% optSolver error add in empty '{}'
tempZono = NN.Z;
leaves = tempZono.getLeaves({});

[ngb, num_leaves] = size(leaves);

% Find a Different way to import faces and vertices
% For larger systems, needing to falsly plot the system is bad
if fromScratch == true
    [LCvertices, LCfacets] = plotHybZono3D(tempZono,{});
    save('vertciesAndFacets.mat',"LCvertices","LCfacets");
else
    load('vertciesAndFacets.mat');
end

[Rdims, Cdims] = size(LCvertices);

%% Calculating Lipschitz constant

% Variable L will contain the largest slope for each facet
L = zeros(num_leaves,1);

% Fix the pipeline, we should store all the vertices as leafs before doing
% any slope calculations should compress the code quite a bit
for leaf = 1:num_leaves
    
    allVerticesForLeaf = [];
    vertexNum = size(rmmissing(LCfacets(leaf,:)),2);
    for i = 1 : vertexNum
        allVerticesForLeaf = [allVerticesForLeaf; LCvertices(LCfacets(leaf,i),:)]; 
    end
    leafRank = rank(allVerticesForLeaf);

%     disp(allVerticesForLeaf)
%     fprintf('Vertices: %d \t Rank: %d\n', vertexNum, leafRank);

    % Quick checks for basic cases
    if (leafRank < 1)
        fprintf('There are no vertices.\n')
    elseif (leafRank == 1)
        fprintf('This is a point.\n');
    elseif (leafRank == 2)
        % Finds the slope for the a plane using the normal vector
        vector = LCvertices(LCfacets(leaf,1),:)-LCvertices(LCfacets(leaf,2),:);

        L(leaf) = abs(vector(3))/ norm(vector(1:2));

        fprintf('This is a line. Slope = %1.4f\n', L(leaf));

    elseif (leafRank == 3)
        vectorOne = LCvertices(LCfacets(leaf,1),:)-LCvertices(LCfacets(leaf,2),:);
        vectorTwo = LCvertices(LCfacets(leaf,1),:)-LCvertices(LCfacets(leaf,3),:);

        normalVector = cross(vectorOne, vectorTwo);

        % Warning gets passed if tolerance is properly set during
        % construction of vertices and facets
        if abs(normalVector(3)) < 1e-6  % "sideways" plane
            warning('Division by 0 in normal vector\n')
        else
            slopeOne = -normalVector(1)/normalVector(3);
            slopeTwo = -normalVector(2)/normalVector(3);

            L(leaf) = norm([slopeOne, slopeTwo]);
        end

        fprintf('This is a plane. Slope = %1.4f\n', L(leaf));
    else    % 3+ vertices
        % Probable Theoretical Pipeline
        % collect all the vertices,
        % find the rank and then re-evaluate the problem again,

        % When more vertices are given issue may arise when selecting
        % points that might be colinear causing adverse effects?

        fprintf('Higher Dimension 3+\n');
    end
end


%% Finding the largest L disregarding outliers and also the outliers
maxL = [];
bigL = [];
outliers = isoutlier(L,"mean");
for i  = 1: length(outliers)
    if(~outliers(i))
        maxL = [maxL, L(i)];
    else
        bigL = [bigL, L(i)];
    end
end
maxAll = max(maxL)
averageL = mean(maxL)

% %% Testing net points,
% load('linear_4_5_4.mat',"net")
% x1 = -1.576 + 0.0010;
% x2 = -1.576 - 0.0010;
% y1 = 5 + 0.0010;
% y2 = 5 - 0.0010;
% points = 21;
% 
% xspan = linspace(x1,x2,points);
% yspan = linspace(y1,y2,points);
% 
% evalPT = [];
% 
% for i = 1:length(xspan)
%     for j = 1: length(yspan)
%         evalPT = [evalPT; xspan(j), yspan(i), net.predict([xspan(j),yspan(i)])];
%     end
% end
% 
% scatter3(evalPT(:,1),evalPT(:,2),evalPT(:,3),'.');
% 
% hold on
% 
% plot3(evalPT(221,1),evalPT(221,2),evalPT(221,3),'o',"Color","k");
