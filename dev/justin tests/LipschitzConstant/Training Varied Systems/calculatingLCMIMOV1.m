% Programmer Name: Justin Chen
% The following program is an attempt to develop a method to calculate the
% lipschtiz constant from trained NN data, MIMO system

% Check version three to prove the calculations for 3D cases. 
clear; close all;
clc;
format short
%% Loading Data
% Loads pretrained data for 2d and 3d cases
% load("S.mat")
% rng(S)

%% Load Neural Net
% load('2in-2out-MIMO_8_4_3',"NN","inDims","outDims","A")
% load('2in-3out-MIMO_8_4_3',"NN","inDims","outDims","A")
load('2in-3out-MIMOV2_8_4_3',"NN","inDims","outDims","A")
% 
leaves = round(NN.Z.getLeaves({}));%tempZono.getLeaves({});%

[ngb, num_leaves] = size(leaves);

waitbarHandle = waitbar(0,['Finding Vertices for the ',num2str(num_leaves),' leaf.']);
optPlot = plotOptions('Display','off','SolverOpts',{});

% Variable L will contain the largest slope for each facet
L = zeros(num_leaves,numel(outDims));
theoryMatrixNorm = zeros(num_leaves,1);
count = 0;

tol = 1e-7;

for leaf = 1 : num_leaves
    for outs = 1: numel(outDims)
        tempZono = NN.projection([inDims outDims(outs)]);
        Zi = conZono(tempZono.Gc,tempZono.c+tempZono.Gb*leaves(:,leaf),tempZono.Ac,tempZono.b-tempZono.Ab*leaves(:,leaf));
    
        try
            [vertices] = findVF(Zi,tol,optPlot.SolverOpts);
        catch
            count = count + 1;
            continue
        end
        
        % Possibly not needed.
        vertices = uniquetol(vertices, tol,"ByRows",true);
        leafRank = rank(vertices);
    
        % Quick checks for basic cases and Calculating LC
        if (leafRank < 1)
            fprintf('There are no vertices.\n')
        elseif (leafRank == 1)
            fprintf('This is a point.\n');
            pause(5)
        elseif (leafRank == 2)
            % Finds the slope for the a plane using the normal vector
            vectorOne = vertices(1,:)-vertices(2,:);
            L(leaf,outs) = abs(vectorOne(3))/ norm(vectorOne(1:2));
            
            calcA(outs,:) = [abs(vectorOne(3))/ norm(vectorOne(1:2)) 0];

            fprintf('This is a line. Slope = %1.4f\n', L(leaf,outs));
            pause(5)
        else
            % Caclulating the vectors based on the found vertices
            for i = 2: leafRank
                vectorOne(i-1,:) = vertices(1,:) - vertices(i,:);
            end
    
            % Find the nullspace using reduce row echelon form
            nullSpace = rref(vectorOne);
            
            calcA(outs,:) = nullSpace(:,leafRank)';
            
            normalVectorRREF = [-nullSpace(:,leafRank); 1];
    
            if abs(normalVectorRREF(leafRank)) < 1e-6  % "sideways" plane
                warning('Division by 0 in normal vector\n')
            else
                for i = 1: leafRank - 1
                    slope(i,:) = -normalVectorRREF(i,:)/normalVectorRREF(leafRank,:);
                end
                L(leaf,outs) = norm(slope);
            end
    
            fprintf('There are %d-points. Slope = %1.4f\n', leafRank, L(leaf,outs));
        end
        waitbar(leaf/num_leaves,waitbarHandle)
    end
    theoryMatrixNorm(leaf) = norm(calcA,2);
end
close(waitbarHandle)

%% Finding the largest L disregarding outliers and also the outliers
avgL = mean(L);

overallL = norm(avgL,2);

worstL = max(vecnorm(L,2,2));

slopeA = norm(vecnorm(A,2,2));

normA = norm(A,2);

theoryNormA = mean(theoryMatrixNorm);


fprintf('\n==========================================================\n');
fprintf('The norm of matrix A is: %1.1f\n', normA);
fprintf('The theoretical norm of matrix A is: %1.1f\n\n', theoryNormA);

fprintf('The slope of matrix A is: %1.1f\n', slopeA);
fprintf('The Average Lipschitz Constant is: %1.1f\n', overallL);


fprintf('The Largest Lipschitz Constant is: %1.1f\n', worstL);