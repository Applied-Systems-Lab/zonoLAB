% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       A dimension-aware and memory-enabled intersection of memZono objects
%   Syntax:
%       [Z] = memorySum(X,Y,sharedDimLabels)
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
%       Shared dimensions undergo an intersection and unshared dimensions 
%       are maintained. Intersections will result in additional constraints 
%       provided and labeled acording to sharedDimLabels.
%       Factors are aligned to preserve memory.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function obj = memorySum(obj1,obj2,sharedDimLabels)
    arguments
        obj1
        obj2
        sharedDimLabels = {};
    end

    %% Memory CartProd
    [Z,keysStruct] = memoryCartProd(obj1,obj2);
    [G_,c_,A_,b_,vset_,keys_] = deal(Z.G_,Z.c_,Z.A_,Z.b_,Z.vset_,Z.keys_);
    lbl = keysStruct.lbl;
    idx = keysStruct.idx;
    
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

    %% Shared Constraints
    if ~isempty(lbl.cs)
        if all([isnumeric(obj1.A_),isnumeric(obj2.A_),isnumeric(obj1.b_),isnumeric(obj2.b_)])
            if all(obj1.A_(idx.cs1,idx.ks1) ~= obj2.A_(idx.cs2,idx.ks2),'all') || all(obj1.b_(idx.cs1) ~= obj2.b_(idx.cs2),'all')
                    error('Shared Constraints are not identical')
            end
        end
        A_ = [A_;
            zeros(length(lbl.cs),length(lbl.k1)), obj1.A_(idx.cs1,idx.ks1), zeros(length(lbl.cs),length(lbl.k2))
        ];
        b_ = [b_;
            obj1.b_(idx.cs1,:)
        ];
        % Labeling
        keys_.cons = [keys_.cons,lbl.cs];
    end

    %% Define memZono
    obj = memZono(G_,c_,A_,b_,vset_,keys_);











    % % Input Conditioning
    % if ~isa(obj1,'memZono') || ~isa(obj2,'memZono')
    %     error('Inputs must both be memZono objects')
    % end
    
    % %% Keys
    % % shared factors
    % [k1,ks,k2] = memZono.getUniqueKeys(obj1.factorKeys,obj2.factorKeys);
    % [idxk1,idxks1,idxks2,idxk2] = memZono.getKeyIndices(obj1.factorKeys,obj2.factorKeys);
    % % shared dims
    % [d1,ds,d2] = memZono.getUniqueKeys(obj1.dimKeys,obj2.dimKeys);
    % [idxd1,idxds1,idxds2,idxd2] = memZono.getKeyIndices(obj1.dimKeys,obj2.dimKeys);
    % % shared cons
    % [c1,cs,c2] = memZono.getUniqueKeys(obj1.conKeys,obj2.conKeys);
    % [idxc1,idxcs1,idxcs2,idxc2] = memZono.getKeyIndices(obj1.conKeys,obj2.conKeys);

    % %% Factor-based Memory Cartisian Product
    % G_ = [
    %     obj1.G_(idxd1,idxk1), obj1.G_(idxd1,idxks1), zeros(length(d1),length(k2));
    %     zeros(length(d2),length(k1)), obj2.G_(idxd2,idxks2), obj2.G_(idxd2,idxk2)
    %     ];
    % c_ = [
    %     obj1.c_(idxd1,:);
    %     obj2.c_(idxd2,:)
    %     ];
    % A_ = [
    %     obj1.A_(idxc1,idxk1), obj1.A_(idxc1,idxks1), zeros(length(c1),length(k2));
    %     zeros(length(c2),length(k1)), obj2.A_(idxc2,idxks2), obj2.A_(idxc2,idxk2);
    %     ];
    % b_ = [
    %     obj1.b_(idxc1,:);
    %     obj2.b_(idxc2,:);
    %     ];

    % % hybrid Zono
    % if obj1.vset_(idxks1) ~= obj2.vset_(idxks2)
    %     error('c/d factors not lining up');
    % end
    % vset_ = [obj1.vset_(idxk1),obj1.vset_(idxks1),obj2.vset_(idxk2)];

    % % Labeling
    % keys_.factors = [k1,ks,k2];
    % keys_.dims = [d1,d2];
    % keys_.cons = [c1,c2];

    % %% Shared Dimensions
    % if ~isempty(ds)
    %     % Intersection matrix definition
    %     R = zeros(length(ds),obj1.n);
    %     for k = 1:length(ds)
    %         j = idxds1(k); 
    %         R(k,j) = 1; 
    %     end
        
    %     % Interesecting Matrices
    %     G_ = [G_;
    %         obj1.G_(idxds1,idxk1), obj1.G_(idxds1,idxks1), zeros(length(ds),length(k2));
    %     ];
    %     c_ = [c_;
    %         obj1.c_(idxds1,:);
    %     ];
    %     A_ = [A_; 
    %         R*obj1.G_(:,idxk1), R*obj1.G_(:,idxks1)-obj2.G_(idxds2,idxks2), -obj2.G_(idxds2,idxk2) 
    %     ];
    %     b_ = [b_;
    %         obj2.c_(idxds2,:) - R*obj1.c_;
    %     ];
    %     % Labels
    %     keys_.dims = [keys_.dims, ds];

    %     % Constraint Labels
    %     if ~isa(sharedDimLabels,'cell')
    %         if ~(isstring(sharedDimLabels)||ischar(sharedDimLabels))
    %             error('Intersection operation will add additional constraints but labels for new constraints are not given.'); 
    %         end
    %         cds{length(ds)} = [];
    %         for k = 1:length(ds)
    %             % TODO: Warning if these new conKey labels already exist in either obj1 or obj2
    %             cds{k} = sprintf('%s_%s_%d',sharedDimLabels,ds{k},k);
    %         end
    %     elseif length(sharedDimLabels) ~= length(ds)
    %         error('Intersection operation will add additional constraints but labels for new constraints are not given (or not enough given).');
    %     else
    %         cds = sharedDimLabels;
    %     end
    %     if ~isempty([obj1.conKeys,obj2.conKeys])
    %         if any(ismember(cds,[obj1.conKeys,obj2.conKeys]))
    %             error('Intersection labels already exists... provide new names')
    %         end
    %     end
    %     keys_.cons = [keys_.cons, cds];
    % end

    % %% Shared Constraints
    % if ~isempty(cs)
    %     if all([isnumeric(obj1.A_),isnumeric(obj2.A_),isnumeric(obj1.b_),isnumeric(obj2.b_)])
    %         if all(obj1.A_(idxcs1,idxks1) ~= obj2.A_(idxcs2,idxks2),'all') || all(obj1.b_(idxcs1) ~= obj2.b_(idxcs2),'all')
    %                 error('Shared Constraints are not identical')
    %         end
    %     end
    %     A_ = [A_;
    %         zeros(length(cs),length(k1)), obj1.A_(idxcs1,idxks1), zeros(length(cs),length(k2))
    %     ];
    %     b_ = [b_;
    %         obj1.b_(idxcs1,:)
    %     ];
    %     % Labeling
    %     keys_.cons = [keys_.cons,cs];
    % end

    % %% Define memZono
    % obj = memZono(G_,c_,A_,b_,vset_,keys_);
    
end
