function obj = cartProd(obj1,obj2,dims1,dims2,options)
    arguments
        obj1 memZono
        obj2 memZono
        dims1 = [];
        dims2 = [];
        options.s1 = '_s1';
        options.s2 = '_s2';
        options.sharedMethod = 'reject';
    end
    if isempty(dims1); dims1 = obj1.dimKeys; end
    if isempty(dims2); dims2 = obj2.dimKeys; end
    [k1,ks,k2] = memZono.getUniqueKeys(dims1,dims2);
    switch options.sharedMethod
        case 'reject'
            if ~isempty(ks)
                error('standard cartProd only works if no dims are in common');
            end
            obj = memoryCartProd(obj1,obj2);
        case 'rename'
            if ~isempty(ks)
                Z1 = obj1.projection(k1); 
                Z2 = obj2.projection(k2);
                Z1s = obj1.projection(ks).relabel(options.s1,'dims');
                Z2s = obj2.projection(ks).relabel(options.s2,'dims');
                % perform cartProd()
                obj = memoryCartProd(Z1s,Z2s);
                if Z1.n>0; obj.memoryCartProd(Z1); end
                if Z2.n>0; obj.memoryCartProd(Z2); end
            else
                obj = memoryCartProd(obj1,obj2);
            end
    end
end