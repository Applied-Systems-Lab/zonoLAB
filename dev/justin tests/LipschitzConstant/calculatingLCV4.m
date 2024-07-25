% Programmer Name: Justin Chen
% The following program is an attempt to develop a method to calculate the
% lipschtiz constant from trained NN data

% Check version three to prove the calculations for 3D cases. 
clear; close all;
clc;

%% Loading last saved data
fromScratch = true;

%% Loading Data
% Loads pretrained data for 2d and 3d cases
load("S.mat")
rng(S)

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
L = zeros(num_leaves,3);

tol = 1e-7;
for leaf = 1 : num_leaves
    Zi{leaf,1} = conZono(tempZono.Gc,tempZono.c+tempZono.Gb*leaves(:,leaf),tempZono.Ac,tempZono.b-tempZono.Ab*leaves(:,leaf));
    [vertices] = findVF(Zi{leaf},tol,optPlot.SolverOpts);
    
    % Possibly not needed.
    vertices = uniquetol(vertices, tol,"ByRows",true);
    leafRank = rank(vertices);

    % Quick checks for basic cases and Calculating LC
    if (leafRank < 1)
        fprintf('There are no vertices.\n')
    elseif (leafRank == 1)
        fprintf('This is a point.\n');
        L(leaf,2) = 1;
    elseif (leafRank == 2)
        % Finds the slope for the a plane using the normal vector
        vectorOne = vertices(1,:)-vertices(2,:);
        L(leaf,1) = abs(vectorOne(3))/ norm(vectorOne(1:2));

        vectorTwo = LCvertices(LCfacets(leaf,1),:)-LCvertices(LCfacets(leaf,2),:);
        L(leaf,3) = abs(vectorTwo(3))/ norm(vectorTwo(1:2));
        
        fprintf('This is a line. Slope = %1.4f\n', L(leaf,1));
        L(leaf,2) = 2;
    else
        % Caclulating the vectors based on the found vertices
        for i = 2: leafRank
            vectorOne(i-1,:) = vertices(1,:) - vertices(i,:);
        end

        % Find the nullspace using reduce row echelon form
        nullSpace = rref(vectorOne);
        normalVectorRREF = [-nullSpace(:,leafRank); 1];

        if abs(normalVectorRREF(leafRank)) < 1e-6  % "sideways" plane
            warning('Division by 0 in normal vector\n')
        else
            for i = 1: leafRank - 1
                slope(i,:) = -normalVectorRREF(i,:)/normalVectorRREF(leafRank,:);
            end
            L(leaf,1) = norm(slope);
        end

        fprintf('There are %d-points. Slope = %1.4f\n', leafRank, L(leaf,1));
        L(leaf,2) = leafRank;
    end
    waitbar(leaf/num_leaves,waitbarHandle)
end
close(waitbarHandle)

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