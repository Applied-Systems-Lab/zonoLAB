% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       A dimension-aware and memory-enabled intersection of memZono objects
%   Syntax:
%       [Z] = memoryIntersection(X,Y,sharedDimLabels)
%   Inputs:
%       X               - memZono in R^n
%       Y               - memZono in R^m
%       sharedDimLabels - either (1) a cell array of new constraint key 
%                                    labels
%                                (2) a string prefix to use for new 
%                                    constraint key labels
%   Outputs:
%       Z - memZono in R^p, where n,m < p <= n+m
%           Shared dimensions are intersected, unshared dimensions kept
%   Notes:
%       - Shared dimensions undergo an intersection and unshared dimensions 
%       are maintained. Intersections will result in additional constraints 
%       provided and labeled acording to sharedDimLabels.
%       - Factors are aligned to preserve memory.
%       - If there are no shared dimensions, memoryIntersection has the same
%       output as cartProd().
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function obj = memoryIntersection(obj1,obj2,sharedDimLabels)
    arguments
        obj1
        obj2
        sharedDimLabels
    end

    %% Memory CartProd
    Z = memoryCartProd(obj1,obj2);
    [G_,c_,A_,b_,vset_,keys_] = Z.exportAllData;
    [lbl,idx] = memZono.getAllKeyIndices(obj1,obj2);
    
    %% Shared Dimensions
    if ~isempty(lbl.ds)
        % Intersection matrix definition
        R = zeros(length(lbl.ds),obj1.n);
        for k = 1:length(lbl.ds)
            j = idx.ds1(k); 
            R(k,j) = 1; 
        end
        
        % Interesecting Matrices
        G_ = [G_;
            obj1.G_(idx.ds1,idx.k1), obj1.G_(idx.ds1,idx.ks1), zeros(length(lbl.ds),length(lbl.k2));
        ];
        c_ = [c_;
            obj1.c_(idx.ds1,:);
        ];
        A_ = [A_; 
            R*obj1.G_(:,idx.k1), R*obj1.G_(:,idx.ks1)-obj2.G_(idx.ds2,idx.ks2), -obj2.G_(idx.ds2,idx.k2) 
        ];
        b_ = [b_;
            obj2.c_(idx.ds2,:) - R*obj1.c_;
        ];
        % Labels
        keys_.dims = [keys_.dims, lbl.ds];

        % Constraint Labels
        if ~isa(sharedDimLabels,'cell')
            if ~(isstring(sharedDimLabels)||ischar(sharedDimLabels))
                error('Intersection operation will add additional constraints but labels for new constraints are not given.'); 
            end
            cds{length(lbl.ds)} = [];
            for k = 1:length(lbl.ds)
                % TODO: Warning if these new conKey labels already exist in either obj1 or obj2
                cds{k} = sprintf('%s_%s_%d',sharedDimLabels,lbl.ds{k},k);
            end
        elseif length(sharedDimLabels) ~= length(lbl.ds)
            error('Intersection operation will add additional constraints but labels for new constraints are not given (or not enough given).');
        else
            cds = sharedDimLabels;
        end
        if ~isempty([obj1.conKeys,obj2.conKeys])
            if any(ismember(cds,[obj1.conKeys,obj2.conKeys]))
                error('Intersection labels already exists... provide new names')
            end
        end
        keys_.cons = [keys_.cons, cds];
    end

    %(shared constraints in memCartProd)
    % %% Shared Constraints
    % if ~isempty(lbl.cs)
    %     if all([isnumeric(obj1.A_),isnumeric(obj2.A_),isnumeric(obj1.b_),isnumeric(obj2.b_)])
    %         if all(obj1.A_(idx.cs1,idx.ks1) ~= obj2.A_(idx.cs2,idx.ks2),'all') || all(obj1.b_(idx.cs1) ~= obj2.b_(idx.cs2),'all')
    %                 error('Shared Constraints are not identical')
    %         end
    %     end
    %     A_ = [A_;
    %         zeros(length(lbl.cs),length(lbl.k1)), obj1.A_(idx.cs1,idx.ks1), zeros(length(lbl.cs),length(lbl.k2))
    %     ];
    %     b_ = [b_;
    %         obj1.b_(idx.cs1,:)
    %     ];
    %     % Labeling
    %     keys_.cons = [keys_.cons,lbl.cs];
    % end

    %% Define memZono
    obj = memZono(G_,c_,A_,b_,vset_,keys_);
end
