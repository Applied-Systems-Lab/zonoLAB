% Projection is defined for internal use - subsref (indexing) is simpilar syntax
function out = projection(obj,dims,options)
    arguments
        obj
        dims
        options.removeExtraFactors = false;
    end
    if ~iscell(dims) % if not already in cell form
        if strcmp(dims,':'), dims = obj.dimKeys;
        else, dims = obj.keysStartsWith(dims).dims;
        end
    end
    [~,idx] = ismember(dims,obj.dimKeys);
    keys_out = obj.keys_; keys_out.dims = dims;
    out = memZono(obj.G_(idx,:),obj.c_(idx,:),obj.A_,obj.b_,obj.vset_,keys_out);   

    if options.removeExtraFactors
        idx = any(vertcat(any(obj.G_),any(obj.A_))); %<= indices where anything is nonzero
        if ~all(idx)
            keys_out.factors = keys_out.factors(idx);
            out = memZono(out.G_(:,idx),out.c_,out.A_(:,idx),out.b_,out.vset_(idx),keys_out);
        end
    end

end