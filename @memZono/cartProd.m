function obj = cartProd(obj1,obj2)
    arguments
        obj1 memZono
        obj2 memZono
    end
    [~,ks,~] = memZono.getUniqueKeys(obj1.dimKeys,obj2.dimKeys);
    if ~isempty(ks)
        error('standard cartProd only works if no dims are in common');
    end
    obj = memoryCartProd(obj1,obj2);
end