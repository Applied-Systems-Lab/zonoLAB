function varargout = subsref(obj, S)
    % Overload of varargout: https://www.mathworks.com/help/matlab/matlab_oop/code-patterns-for-subsref-and-subsasgn-methods.html
    if length(S) > 1
        switch S(1).type
            case '{}'
                switch S(2).type
                    case {'()','.'}
                        [varargout{1:nargout}] = subsref(subsref(obj,S(1)),S(2:end));
                        return;
                end
            case '()'
                if strcmp(S(2).type,'.')
                    [varargout{1:nargout}] = subsref(subsref(obj,S(1)),S(2:end));
                    return;
                end
            case '.'
                % for i = 2:length(S) %<== allow multiple function calls/properties (should be automatic so commented out)
                %     if strcmp(S(i).type,'.')
                %         [varargout{1:nargout}] = subsref(subsref(obj,S(1:(i-1))),S(i:end));
                %         return;
                %     end
                % end
        end
    end
    switch S(1).type
        case '{}'
            [varargout{1:nargout}] = builtin('subsref', obj, S);
        case '.'
            % if ismember(S(1).subs,{'boundingBox','plot','bounds','lb','ub'}) %<=== explicitly calls these specific methods
            %     [varargout{1:nargout}] = builtin('subsref', Z(obj,obj.dimKeys), S); %<== calls w/ assumption that dims are originally indexed
            %     for i = 1:numel(varargout)
            %         varargout{i} = addDimKeys(varargout{i},obj,S(1).subs);
            %     end
            % else
                [varargout{1:nargout}] = builtin('subsref', obj, S);
            % end
        case '()'
            switch numel(S(1).subs)
                case 1
                    i = S(1).subs{1}; j = ':'; k = ':';
                case 2
                    i = S(1).subs{1}; j = S(1).subs{2}; k = ':';
                case 3
                    i = S(1).subs{1}; j = S(1).subs{2}; k = S(1).subs{3};
                otherwise
                    error('indexing not specificed')
            end

            i = getKeyIndices(i,obj.keys_.dims);
            j = getKeyIndices(j,obj.keys_.factors);
            k = getKeyIndices(k,obj.keys_.cons);

            if any(~[i,j,k])
                if any(~i)
                    error("dims specified for indexing don't all exist");
                else
                    error('indexing refering to non-existant factors/dims');
                end
            end

            G_ = obj.G_(i,j);
            c_ = obj.c_(i,:);
            A_ = obj.A_(k,j);
            b_ = obj.b_(k,:);
            vset_ = obj.vset_(1,j);

            keys_.dims = obj.keys_.dims(i);
            keys_.factors = obj.keys_.factors(j);
            keys_.cons = obj.keys_.cons(k);

            varargout{1} = memZono(G_,c_,A_,b_,vset_,keys_);



            % i = getKeyIndices(i,obj.dimKeys);
            % j = getKeyIndices(j,obj.factorKeys);
            % k = getKeyIndices(k,obj.conKeys);
            % if ischar(i)
            %     keys_.dims = obj.dimKeys;
            % else
            %     keys_.dims = obj.dimKeys(i);
            % end
            % if ischar(j)
            %     keys_.factors = obj.factorKeys;
            % else
            %     keys_.factors = obj.factorKeys(j);
            % end
            % if ischar(k)
            %     keys_.cons = obj.conKeys;
            % else
            %     keys_.cons = obj.conKeys(k);
            % end

            % varargout{1} = memZono(G_,c_,A_,b_,vset_,keys_);
    end
end


function idx = getKeyIndices(in,keys)
    if isnumeric(in)
        idx = in;
    else
        if iscell(in)
            [~,idx] = ismember(in,keys);
        else
            if strcmp(in,':')
                idx = 1:numel(keys);
            else
                [~,idx] = ismember(in,keys);
                if idx == 0 
                    idx = find(cellfun(@(lbl) startsWith(lbl,in), keys)); 
                end
            end
        end
    end
end

% function out = addDimKeys(in,obj,lbl)
%     if isa(in,'abstractZono')
%         sz = [obj.n,obj.nG,obj.nC];
%         sz_ = [in.n,in.nG,max([0,in.nC])];
%         if sz_ == sz %<== maintains same keys
%             out = memZono(in,obj.keys_);
%         elseif sz_ == [sz(1),sz(1),0]
%             obj.keys_.factors = cellfun(@(key) strcat(key,'_',lbl),obj.keys_.dims,UniformOutput=false); %<== factors same as dims for output
%             obj.keys_.cons = [];
%             out = memZono(in,obj.keys_);
%         else
%             out = memZono(in,obj.keys_.dims,'noMem'); %<== factors/constraint results lost
%         end
%     elseif size(in,1) == obj.n
%         out = array2table(in,RowNames=obj.keys_.dims);
%     else
%         warning('dimension maintance lost with method result')
%     end
% end