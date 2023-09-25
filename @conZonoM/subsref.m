function B = subsref(A, S)
    switch S(1).type
        case {'.','{}'}
            B = builtin('subsref', A, S);
        case '()'
            switch numel(S(1).subs)
                case 1
                    i = S(1).subs{1}; j = ':'; k = ':';
                case 2
                    i = S(1).subs{1}; j = S(1).subs{2}; k = ':';
                case 3
                    i = S(1).subs{1}; j = S(1).subs{2}; k = ':';
                otherwise
                    error('indexing not specificed')
            end
            % ensure numeric
            if isnumeric(i); i = A.dimKeys(i); end
            if isnumeric(j); j = A.factorKeys(j); end
            if isnumeric(k); k = A.conKeys(k); end

            % select parts to keep
            B = A.copy();
            B.c_dict = A.c_dict(i,1);
            B.G_dict = A.G_dict(i,j);
            B.A_dict = A.A_dict(k,j);
            B.b_dict = A.b_dict(k,1);
    end
end