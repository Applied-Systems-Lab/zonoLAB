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
function [v] = findVF2(obj,tol,optSolver)

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
            if (~isempty(foundVerts))
                temp = null([foundVerts;zeros(2,3)]);
                searchDir(i,:) = temp(:,1)';
            else
                searchDir(i,:) = dir(indx, :);
            end
            searchDir(i,:) = normalize(searchDir(i,:) ,'norm');
            %[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,optSolver);
            x = findVertex(searchDir(i,:),obj.G,Aeq,beq,lb,ub,optSolver);
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


%% Local Functions
function [x] = findVertex(dir,G,Aeq,beq,lb,ub,optSolver)
    [x,~,~] = solveLP(dir*G,[],[],Aeq,beq,lb,ub,optSolver);
    if isnan(x)
        error('PlotError:VertexNotFound','Could not find a solution for a vertex while plotting')
    end
end