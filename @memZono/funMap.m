function out = funMap(obj1,obj2,inDims,outDims,lbl,funInDims,funOutDims)
    arguments
        obj1 memZono %<== input
        obj2 memZono %<== funciton
        inDims
        outDims
        lbl
        funInDims = {}; %<== same as inDims if not provided
        funOutDims = {}; %<== same as outDims if not provided
    end

    % Keys setup
    if ~iscell(inDims)
        inDims = obj1.keysStartsWith(inDims).dims;
    end
    if ~iscell(outDims)
        outDims = obj1.keysStartsWith(outDims).dims;
    end
    if isempty(funInDims)
        funInDims = inDims;
    elseif ~iscell(funInDims)
        funInDims = obj2.keysStartsWith(funInDims).dims;
    end
    if isempty(funInDims)
        funInDims = inDims;
    else
        if ~iscell(funInDims)
            funInDims = obj2.keysStartsWith(funInDims).dims;
        end
        if numel(funInDims) ~= numel(inDims)
            error("funInDims does not have same number of dimensions as inDims");
        end
    end
    if isempty(funOutDims)
        funOutDims = outDims;
    else
        if ~iscell(funOutDims)
            funOutDims = obj2.keysStartsWith(funOutDims).dims;
        end
        if numel(funOutDims) ~= numel(outDims)
            error("funOutDims does not have same number of dimensions as inDims");
        end
    end

    % Call intersection then projection
    out = projection(and(...
        projection(obj1,inDims),...
        relabelDims(obj2,[funInDims,funOutDims],[inDims,outDims]), ... %<== relabel if different fun dims
    lbl),outDims);
    

end