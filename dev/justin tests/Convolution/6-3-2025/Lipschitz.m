% Programmer Name: Justin Chen
% Final Version of Lipschitz Calculating
% Last Updated: 6-4-2025

clear; close all;
clc;

%% Loading Data
% load('6_2sincos2x_20_10_10.mat')
load('lipschitzminstV2.mat')
NN = stackedNN;

% Current zonoLab version does not have the resolved file, to get around
% optSolver error add in empty '{}'
tempZono = NN.Z(NN.dimKeys);
leaves = tempZono.getLeaves({});

[ngb, num_leaves] = size(leaves);

waitbarHandle = waitbar(0,['Finding Vertices for the ',num2str(num_leaves),' leaf.']);
optPlot = plotOptions('Display','off');

% Variable L will contain the largest slope for each facet
L = zeros(num_leaves,2);

for leaf = 1 : num_leaves
    Zi = conZono(tempZono.Gc,tempZono.c+tempZono.Gb*leaves(:,leaf),tempZono.Ac,tempZono.b-tempZono.Ab*leaves(:,leaf));
    [vertcies] = findVF(Zi,1e-7,optPlot.SolverOpts);
    
    % Possibly not needed.
    vertcies = uniquetol(vertcies, 1e-10,"ByRows",true);

    leafRank = rank(vertcies);

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

        fprintf('This is a line. Slope = %1.4f\n', L(leaf,1));
        L(leaf,2) = 2;
    elseif (leafRank == 3)
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
                subD(i,j) = det(tempVertices);
            end
        end

        fprintf('Higher Dimension 3+ at Leaf: %d\n', leaf);
    end
    waitbar(leaf/num_leaves,waitbarHandle)
end
close(waitbarHandle)