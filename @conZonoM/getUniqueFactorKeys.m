function [k1,kshared,k2] = getUniqueFactorKeys(obj1,obj2)
    if isempty(obj1.factorKeys) || isempty(obj2.factorKeys)
        kshared = {};
        k1 = obj1.factorKeys;
        k2 = obj2.factorKeys;
    else
        kshared = intersect(obj1.factorKeys,obj2.factorKeys);
        k1 = setdiff(obj1.factorKeys,kshared);
        k2 = setdiff(obj2.factorKeys,kshared);
    end
end