function out = minkowskiSum(obj1,obj2)
% Minkowski Sum between two Conzono Objects
arguments
    obj1     conZonoM
    obj2     conZonoM
end

% Geting Unique Factors
[k1,kshared,k2] = getUniqueFactorKeys(obj1,obj2);
out = conZonoM;
out.factorKeys = [k1,kshared,k2];

if obj1.n ~= obj2.n; error('dimension error'); end
out.dimKeys = ensureKeysMatch(obj1.dimKeys,obj2.dimKeys);


out.conKeys = [obj1.conKeys, obj2.conKeys];

% Performing Minkowski Sum operaions
out.c = obj1.c + obj2.c;
out.G_dict = [
    obj1.G_dict(:,k1), ...
    obj1.G_dict(:,kshared) + obj2.G_dict(:,kshared),...
    obj2.G_dict(:,k2)
    ];
out.A_dict = [
    obj1.A_dict(:,k1), obj1.A_dict(:,kshared), zeroTable(obj1.nC,k2);
    zeroTable(obj2.nC,k1),obj2.A_dict(:,kshared), obj2.A_dict(:,k2)
    ];
out.b = [obj1.b; obj2.b];

end

function out = ensureKeysMatch(k1,k2)
    if any([isempty(k1),isempty(k2)])
        if ~isempty(k1); out=k1; 
        elseif ~isempty(k2); out=k2; 
        else; out = [];
        end
    elseif any(cell2mat(cellfun(@(x,y) x~=y, k1, k2, 'UniformOutput', false)))
        warning('dims not lining up')
    else
        out = k1; % since k1 = k2
    end
end