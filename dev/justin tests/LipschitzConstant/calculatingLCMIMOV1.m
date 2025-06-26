% Programmer Name: Justin Chen
% The following program is an attempt to develop a method to calculate the
% lipschtiz constant from trained NN data, MIMO system

% Check version three to prove the calculations for 3D cases. 
clear; close all;
clc;
format short
%% Load Neural Net
% inDims are the dim keys for the input dimensions
% outDims are the dim keys for the out dimensions
% Replace name below
load('mimoNNforClass1.mat',"NN","inDims","outDims")
% load('mimoNN.mat',"NN","inDims","outDims")
% 
leaves = round(NN.Z.getLeaves({}));%tempZono.getLeaves({});%

[ngb, num_leaves] = size(leaves);

waitbarHandle = waitbar(0,['Finding Vertices for the ',num2str(num_leaves),' leaf.']);
optPlot = plotOptions('Display','off','SolverOpts',{});

% Variable L will contain the largest slope for each facet
L = zeros(num_leaves,numel(outDims));
averageA = 0;
count = 0;

tol = 1e-7;

for leaf = 1 : num_leaves
    for outs = 1: numel(outDims)
        tempZono = NN.projection([inDims outDims(outs)]);
        Zi{leaf,1} = conZono(tempZono.Gc,tempZono.c+tempZono.Gb*leaves(:,leaf),tempZono.Ac,tempZono.b-tempZono.Ab*leaves(:,leaf));
        [LCvertices, LCfacets] = plot(Zi{leaf,1},optPlot); 

        % If an empty leaf is found, skip this iteration
        if(Zi{leaf,1}.checkEmpty == 1)
            area(leaf) = 0;
            tempSlope(leaf,:) = zeros(numel(inDims)+numel(outDims),1); % place holder for empty leaves not sure size
            calcA = 0;
            continue
        end

        try
            [vertices] = findVF(Zi{leaf,1},tol,optPlot.SolverOpts);
        catch
            count = count + 1;
            continue
        end
        
        % Possibly not needed.
        vertices = uniquetol(vertices, tol,"ByRows",true);
        area(leaf) = polyarea(LCvertices(:,1),LCvertices(:,2));% for 2d only
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
                tempSlope(leaf,:) = slope;
                L(leaf,outs) = norm(slope);
            end
            fprintf('There are %d-points. Slope = %1.4f\n', leafRank, L(leaf,outs));
        end
        waitbar(leaf/num_leaves,waitbarHandle)
    end
    adjustedL(leaf,:) = (L(leaf,:)*area(leaf))/4;
    theoryMatrixNorm{leaf,1} = {calcA*area(leaf)/4};% 4 is input domain for specific case
end
close(waitbarHandle)

%% Structure of L

% the first m elements of out are the slopes calculated from the outdims. 

%% Finding the largest L disregarding outliers and also the outliers
avgL = mean(L);

overallL = norm(avgL,2);

worstL = max(vecnorm(L,2,2));



fprintf('\n==========================================================\n');

fprintf('The Average Lipschitz Constant is: %1.1f\n', overallL);

fprintf('The Largest Lipschitz Constant is: %1.1f\n', worstL);

%% Local Functions

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Return vertices and faces for a constrained zonotope in higher
%       dimensions "4D"
%   Syntax:
%       [v,f] = findVertex(Z,optSolver)
%   Inputs:
%       Z - higher dimension constrained zonotope in CG-Rep (conZono object)
%       optSolver - solver options needed for linear propgram
%   Outputs:
%       v - nV x m matrix, each row denoting the x (first column) and y (second column) positions
%                          of the nV vertices
%       f - 1 x nV vector, indicating a single face containing all nV vertices 
%   Notes:
%       Not intended to be called directly by user.
%       
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [v] = findVF(obj,tol,optSolver)

% Problem data for linear program (LP)
Aeq = sparse(obj.A);
beq = [obj.b];
lb = -ones(obj.nG,1);
ub =  ones(obj.nG,1);

%     % Loop through each basis vector to find a vertex in that direction
dir = [  eye(obj.n);
        -eye(obj.n)];
rankVerts = 1; indx = 1;
for j = 0 : 1000
    % Redefine dir based obj.n each loop for each direction
    for i = rankVerts : obj.n
        % Find first vertex
        if indx >= length(dir)
            searchDir = 2*rand(1,obj.n) - 1;
        else
            searchDir = dir(indx, :);
        end
        searchDir = normalize(searchDir ,'norm');
        %[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,optSolver);
        x = findVertex(searchDir,obj.G,Aeq,beq,lb,ub,optSolver);
        foundVerts(i,:) = [obj.G*x + obj.c]';
        indx = indx + 1;
    end

    if rank(foundVerts,tol) == obj.n
        break;
    else
        foundVerts = uniquetol(foundVerts, tol,"ByRows",true);
        rankVerts = rank(foundVerts,tol) + 1;
    end
end

%     if rank(foundVerts) ~= obj.n
%         warning('rank of %d not found', obj.n);
%     end

v = foundVerts;
return;
end

function [x] = findVertex(dir,G,Aeq,beq,lb,ub,optSolver)
[x,~,~] = solveLP(dir*G,[],[],Aeq,beq,lb,ub,optSolver);
if isnan(x)
    error('PlotError:VertexNotFound','Could not find a solution for a vertex while plotting')
end
end