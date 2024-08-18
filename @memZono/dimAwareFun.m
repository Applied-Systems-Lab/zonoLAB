function varargout = dimAwareFun(obj,fun,dimIn,dimOut,lbl,options)
    arguments
        obj memZono
        fun function_handle
        dimIn = [];
        dimOut = [];
        lbl = [];
        options.flag = []; %<== no flags used yet
    end

    if isempty(dimIn); dimIn = obj.dimKeys; end
    if isempty(dimOut); dimOut = obj.dimKeys; end
    if isempty(lbl); lbl = func2str(fun); end
    [varargout{1:nargout}] = feval(fun,obj.Z(dimIn)); %<== eval as baseZono version
    for i = 1:numel(varargout)
        varargout{i} = addDimKeys(varargout{i},dimOut,lbl);
    end
end

function out = addDimKeys(in,dimOut,lbl)
    if isa(in,'abstractZono') % base-zono out
        out = memZono(in,dimOut,lbl);
        if out.nG == out.n
            out.factorKeys = cellfun(...
                @(key)strcat(key,'_',lbl),out.dimKeys,UniformOutput=false);
        end
        if out.nC == out.n
            out.conKeys = cellfun(...
                @(key)strcat(key,'_',lbl),out.conKeys,UniformOutput=false);
        end
    else
    % elseif ~isscalar(in) % non-sclar array
        outDims = memZono.keysCheck(dimOut,size(in,1));
        out = array2table(in,RowNames=outDims);
    % elseif isscalar(in)
    %     warning("scalar output ignores dimesnions")
    % else
    %     error('dimension maintance lost with method result')
    end
end