function varargout = dimAwareFun(fun,varargin,options)
    arguments
        fun function_handle
    end
    arguments (Repeating)
        varargin memZono
    end
    arguments
        options.dimIn = [];
        options.dimOut = [];
        options.lbl = [];
        options.flag = []; %<== no flags used yet
    end

    % dims based on original input:
    if isscalar(varargin)
        obj = varargin{1};
        % Get defaults
        if isempty(options.dimIn); options.dimIn = obj.dimKeys; end
        if isempty(options.dimOut); options.dimOut = obj.dimKeys; end
        if isempty(options.lbl); options.lbl = func2str(fun); end
        % Eval Fun
        [Z_{1:nargout}] = feval(fun,obj.Z(options.dimIn)); %<== eval as baseZono version
        [varargout{1:nargout}] = cellfun(@(Z) addDimKeys(Z,options.dimOut,options.lbl),Z_);
        % for i = 1:numel(varargout)
        %     varargout{i} = addDimKeys(Z_{i},dimOut,lbl);
        % end
    else
        if any([isempty(options.dimIn),isempty(options.dimOut),isempty(options.lbl)])
            error('need to supply all inDims/outDims/lbl for multiple input functions')
        end
        % pre-process
        Z_ = cellfun(...
            @(in,dimIn)in.Z(dimIn),varargin,options.dimIn,...
            UniformOutput=false);
        % eval fun
        [out_{1:nargout}] = feval(fun,Z_{:}); %<== eval as baseZono version
        % post-process
        varargout = cellfun(...
            @(out) addDimKeys(out,options.dimOut,options.lbl),out_,...
            UniformOutput=false);
    end
end


% %% original
% function varargout = dimAwareFun(fun,obj,dimIn,dimOut,lbl,options)
%     arguments
%         fun function_handle
%         obj memZono
%         dimIn = [];
%         dimOut = [];
%         lbl = [];
%         options.flag = []; %<== no flags used yet
%     end

%     if isempty(dimIn); dimIn = obj.dimKeys; end
%     if isempty(dimOut); dimOut = obj.dimKeys; end
%     if isempty(lbl); lbl = func2str(fun); end
%     [varargout{1:nargout}] = feval(fun,obj.Z(dimIn)); %<== eval as baseZono version
%     for i = 1:numel(varargout)
%         varargout{i} = addDimKeys(varargout{i},dimOut,lbl);
%     end
% end

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
        outDims = memZono.keysCheck(dimOut,size(in,1));
        out = array2table(in,RowNames=outDims);
    end
end





