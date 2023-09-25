function [k1,kshared,k2] = getUniqueDimKeys(obj1,obj2)
    if isempty(obj1.dimKeys) || isempty(obj2.dimKeys)
        kshared = {};
        k1 = obj1.dimKeys;
        k2 = obj2.dimKeys;
    else
        kshared = intersect(obj1.dimKeys,obj2.dimKeys);
        k1 = setdiff(obj1.dimKeys,kshared);
        k2 = setdiff(obj2.dimKeys,kshared);
    end
end