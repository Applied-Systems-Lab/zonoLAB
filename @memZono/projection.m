% Projection is defined for internal use - subsref (indexing) is simpilar syntax
function out = projection(obj,dims)
    if ~iscell(dims) % if not already in cell form
        if strcmp(dims,':'), dims = obj.dimKeys;
        else, dims = obj.keysStartsWith(dims).dims;
        end
    end
    [~,idx] = ismember(dims,obj.dimKeys);
    keys_out = obj.keys_; keys_out.dims = dims;
    out = memZono(obj.G_(idx,:),obj.c_(idx,:),obj.A_,obj.b_,obj.vset_,keys_out);
end