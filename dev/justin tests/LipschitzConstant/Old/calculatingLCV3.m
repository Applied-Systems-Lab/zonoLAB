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

waitbarHandle = waitbar(0,['Finding Vertices for the ',num2str(num_leaves),' leaf.']);
optPlot = plotOptions('Display','off','SolverOpts',{});

% Variable L will contain the largest slope for each facet
L = zeros(num_leaves,2);

for leaf = 1 : num_leaves
    Zi = conZono(tempZono.Gc,tempZono.c+tempZono.Gb*leaves(:,leaf),tempZono.Ac,tempZono.b-tempZono.Ab*leaves(:,leaf));
    [vertcies] = findVF(Zi,optPlot.SolverOpts);
    
    % Possibly not needed.
    vertcies = uniquetol(vertcies, 1e-10,"ByRows",true);

    leafRank = rank(vertcies);
    if leaf == 12 || leaf == 14
        continue;
    end
    % Quick checks for basic cases and Calculating LC
    if (leafRank < 1)
        fprintf('There are no vertices.\n')
    elseif (leafRank == 1)
        fprintf('This is a point.\n');
        L(leaf,2) = 1;
    elseif (leafRank == 2)
        % Finds the slope for the a plane using the normal vector
        vector = vertcies(1,:)-vertcies(2,:);

        L(leaf,1) = abs(vector(3))/ norm(vector(1:2));
        Ln(leaf,1) = L(leaf,1);

        fprintf('This is a line. Slope = %1.4f\n', L(leaf,1));
        L(leaf,2) = 2;
    elseif (leafRank == 3)
        % Caclulating the vectors based on the found vertices
        for i = 2: leafRank
            vector(i-1,:) = vertcies(1,:) - vertcies(i,:);
        end

        % Find the nullspace using reduce row echelon form
        nullSpace = rref(vector);
        normalVectorRREF = [-nullSpace(:,3); 1];

        if abs(normalVectorRREF(3)) < 1e-6  % "sideways" plane
            warning('Division by 0 in normal vector\n')
        else
            for i = 1: leafRank - 1
                slope(i,:) = -normalVectorRREF(i,:)/normalVectorRREF(leafRank,:);
            end
            Ln(leaf,1) = norm(slope);
        end


        vectorOne = vertcies(1,:)-vertcies(2,:);
        vectorTwo = vertcies(1,:)-vertcies(3,:);

        normalVector = cross(vectorOne, vectorTwo);

        % Warning gets passed if tolerance is properly set during
        % construction of vertices and facets
        if abs(normalVector(3)) < 1e-6  % "sideways" plane
            warning('Division by 0 in normal vector\n')
        else
            slopeOne = -normalVector(1)/normalVector(3);
            slopeTwo = -normalVector(2)/normalVector(3);

            L(leaf,1) = norm([slopeOne, slopeTwo]);
        end

        fprintf('This is a plane. Slope = %1.4f\n', L(leaf,1));
        L(leaf,2) = 3;
    else    % 3+ vertices
        subD = zeros(length(vertcies));
        for i = 1 : length(vertcies)
            for j = 1 : length(vertcies)
                tempVertices = vertcies([1:i-1 i+1:length(vertcies)], [1:j-1 j+1:length(vertcies)]);
%                 tempVertices = tempVertices(:, [1:j-1 j+1:length(vertcies)]);
                subD(i,j) = det(tempVertices);
            end
        end

%         Li = norm(m);
%         L(leaf,1) = Li;
%         L(leaf,2) = length(vertices);
        fprintf('Higher Dimension 3+ at Leaf: %d\n', leaf);
    end
    waitbar(leaf/num_leaves,waitbarHandle)
end
close(waitbarHandle)

for i = 1: length(Ln)
    if(L(i,2)==3)
        if(abs(Ln(i) - L(i,1)) > 1e-8)
            warning('L and Ln are not the same');
        end
    end
end

% %% Old LC with Hyb
% if fromScratch == true
%     [LCvertices, LCfacets] = plotHybZono3D(tempZono,{});
%     save('vertciesAndFacets.mat',"LCvertices","LCfacets");
% else
%     load('vertciesAndFacets.mat');
% end
% 
% [Rdims, Cdims] = size(LCvertices);
% 
% %% Calculating Lipschitz constant
% 
% % Variable L will contain the largest slope for each facet
% Lo = zeros(num_leaves,2);
% 
% % Fix the pipeline, we should store all the vertices as leafs before doing
% % any slope calculations should compress the code quite a bit
% for leaf = 1:num_leaves
%     
%     allVerticesForLeaf = [];
%     vertexNum = size(rmmissing(LCfacets(leaf,:)),2);
%     for i = 1 : vertexNum
%         allVerticesForLeaf = [allVerticesForLeaf; LCvertices(LCfacets(leaf,i),:)]; 
%     end
%     leafRank = rank(allVerticesForLeaf);
% 
% %     disp(allVerticesForLeaf)
% %     fprintf('Vertices: %d \t Rank: %d\n', vertexNum, leafRank);
% 
%     % Quick checks for basic cases
%     if (leafRank < 1)
%         fprintf('There are no vertices.\n')
%     elseif (leafRank == 1)
%         fprintf('This is a point.\n');
%         Lo(leaf,2) = 1;
%     elseif (leafRank == 2)
%         % Finds the slope for the a plane using the normal vector
%         vector = LCvertices(LCfacets(leaf,1),:)-LCvertices(LCfacets(leaf,2),:);
% 
%         Lo(leaf,1) = abs(vector(3))/ norm(vector(1:2));
% 
%         fprintf('This is a line. Slope = %1.4f\n', Lo(leaf,1));
%         Lo(leaf,2) = 2;
% 
%     elseif (leafRank == 3)
%         vectorOne = LCvertices(LCfacets(leaf,1),:)-LCvertices(LCfacets(leaf,2),:);
%         vectorTwo = LCvertices(LCfacets(leaf,1),:)-LCvertices(LCfacets(leaf,3),:);
% 
%         normalVector = cross(vectorOne, vectorTwo);
% 
%         % Warning gets passed if tolerance is properly set during
%         % construction of vertices and facets
%         if abs(normalVector(3)) < 1e-6  % "sideways" plane
%             warning('Division by 0 in normal vector\n')
%         else
%             slopeOne = -normalVector(1)/normalVector(3);
%             slopeTwo = -normalVector(2)/normalVector(3);
% 
%             Lo(leaf,1) = norm([slopeOne, slopeTwo]);
%         end
% 
%         fprintf('This is a plane. Slope = %1.4f\n', Lo(leaf,1));
%         Lo(leaf,2) = 3;
%     else    % 3+ vertices
%         % Probable Theoretical Pipeline
%         % collect all the vertices,
%         % find the rank and then re-evaluate the problem again,
% 
%         % When more vertices are given issue may arise when selecting
%         % points that might be colinear causing adverse effects?
% 
%         fprintf('Higher Dimension 3+\n');
%     end
% end
% 
% for i = 1: length(Lo)
%     if(Lo(i) - L(i)) > 1e-8 && (Lo(i,2) == L(i,2))
%         fprintf("Different slope at leaf: %d\n",i);
%     end
% end
% 
% for i = 1: length(Lo)
%     if(Lo(i,2) ~= L(i,2))
%         fprintf("Different classification at leaf: %d\n",i);
%     end
% end
% %% Finding the largest L disregarding outliers and also the outliers
% maxL = [];
% bigL = [];
% outliers = isoutlier(L,"mean");
% for i  = 1: length(outliers)
%     if(~outliers(i))
%         maxL = [maxL, L(i)];
%     else
%         bigL = [bigL, L(i)];
%     end
% end
% maxAll = max(maxL)
% averageL = mean(maxL)