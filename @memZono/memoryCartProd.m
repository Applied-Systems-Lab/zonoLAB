% Returns 
function [obj,keysStruct] = memoryCartProd(obj1,obj2)
    arguments
        obj1 memZono
        obj2 memZono
    end

    %% Keys
    % shared dims
    [lbl.d1,lbl.ds,lbl.d2] = memZono.getUniqueKeys(obj1.dimKeys,obj2.dimKeys);
    [idx.d1,idx.ds1,idx.ds2,idx.d2] = memZono.getKeyIndices(obj1.dimKeys,obj2.dimKeys);
    % shared factors
    [lbl.k1,lbl.ks,lbl.k2] = memZono.getUniqueKeys(obj1.factorKeys,obj2.factorKeys);
    [idx.k1,idx.ks1,idx.ks2,idx.k2] = memZono.getKeyIndices(obj1.factorKeys,obj2.factorKeys);
    % shared cons
    [lbl.c1,lbl.cs,lbl.c2] = memZono.getUniqueKeys(obj1.conKeys,obj2.conKeys);
    [idx.c1,idx.cs1,idx.cs2,idx.c2] = memZono.getKeyIndices(obj1.conKeys,obj2.conKeys);

    %% Factor-based Memory Cartisian Product
    G_ = [
        obj1.G_(idx.d1,idx.k1), obj1.G_(idx.d1,idx.ks1), zeros(length(lbl.d1),length(lbl.k2));
        zeros(length(lbl.d2),length(lbl.k1)), obj2.G_(idx.d2,idx.ks2), obj2.G_(idx.d2,idx.k2)
        ];
    c_ = [
        obj1.c_(idx.d1,:);
        obj2.c_(idx.d2,:)
        ];
    A_ = [
        obj1.A_(idx.c1,idx.k1), obj1.A_(idx.c1,idx.ks1), zeros(length(lbl.c1),length(lbl.k2));
        zeros(length(lbl.c2),length(lbl.k1)), obj2.A_(idx.c2,idx.ks2), obj2.A_(idx.c2,idx.k2);
        ];
    b_ = [
        obj1.b_(idx.c1,:);
        obj2.b_(idx.c2,:);
        ];

    % hybrid Zono
    if obj1.vset_(idx.ks1) ~= obj2.vset_(idx.ks2)
        error('c/d factors not lining up');
    end
    vset_ = [obj1.vset_(idx.k1),obj1.vset_(idx.ks1),obj2.vset_(idx.k2)];

    % Labeling
    keys_.factors = [lbl.k1,lbl.ks,lbl.k2];
    keys_.dims = [lbl.d1,lbl.d2];
    keys_.cons = [lbl.c1,lbl.c2];

    %% Define memZono
    obj = memZono(G_,c_,A_,b_,vset_,keys_);

    %% out struct    
    keysStruct.lbl = lbl;
    keysStruct.idx = idx;
    % keysStruct.lbl = cell2struct(...
    %     {lbl.d1,lbl.ds,lbl.d2,...
    %      lbl.k1,lbl.ks,lbl.k2,...
    %      lbl.c1,lbl.cs,lbl.c2},...
    %     {'d1','ds','d2',...
    %      'k1','ks','k2',...
    %      'c1','cs','c2'},2);
    % keysStruct.idx = cell2struct(...
    %     {idx.d1,idx.ds1,idx.ds2,idx.d2,...
    %      idx.k1,idx.ks1,idx.ks2,idx.k2,...
    %      idx.c1,idx.cs1,idx.cs2,idx.c2 ...
    %     },...
    %     {'d1','ds1','ds2','d2',...
    %      'k1','ks1','ks2','k2',...
    %      'c1','cs1','cs2','c2'...
    %     },2);
        % {idx.d1,idx.ds1,idx.ds2,idx.d2,...
        %  idx.k1,idx.ks1,idx.ks2,idx.k2,...
        %  idx.c1,idx.cs1,idx.cs2,idx.c2 ...
        % },...
        % {'idx.d1','idx.ds1','idx.ds2','idx.d2',...
        %  'idx.k1','idx.ks1','idx.ks2','idx.k2',...
        %  'idx.c1','idx.cs1','idx.cs2','idx.c2'...
        % },2);
end