classdef conZonoM < conZono
    % A class for storing constrained zonotopes alongside labels for
    % generators to utilize memory

    % Constrained zonotope of the form
    % X = { c + G \xi | ||\xi||_inf <= 1, A \xi = b }
    % 
    % factorKeys = labeled collumns of G and A
    % dimKeys = labeled rows of c and G
    % conKeys = labeled rows of A and b
    
    % All ConZono properties inherited... additonal properties:
    properties
        factorKeys % keys associated with the factors (collumns of G and A)
        dimKeys % keys associated with x dimensions (rows of c and G)
        conKeys % keys associated with constraints (rows of A and b)
    end

    properties (Dependent, Hidden)
        c_dict % Dict of c
        G_dict % Dict of G
        A_dict % Dict of A
        b_dict % Dict of b
        Z % conZono Object
        keys % The struct of all keys
    end

    methods
        function obj = conZonoM(varargin)
            % conZonoM Constructor

            % Define conZono as a superClass
            obj = obj@conZono();            
            if nargin == 1 || nargin == 2
                % in = varargin{1};
                % obj = conZono(in.G,in.c,in.A,in.b);
                % obj.G = varargin{1}.G;
                % obj.c = varargin{1}.c;
                % obj.A = varargin{1}.A;
                % obj.b = varargin{1}.b;
                obj.Z = varargin{1};
            elseif nargin > 2
                obj.Z = conZono(varargin{1:end-1});
            end

            % Set keys
            if nargin >= 2
                if isa(varargin(end),"struct")
                    obj.keys = varargin(end);
                else % only factors
                    obj.factorKeys = varargin(end);
                    obj.dimKeys = [];
                    obj.conKeys = [];
                end
            end
        end

        %% Getter functions
        function value = get.c_dict(obj)
            % Calculate the value of c_dict based on c and keys
            value = array2table(obj.c);
            if ~isempty(obj.dimKeys)
                value.Properties.RowNames = obj.dimKeys;
            end
        end

        function value = get.G_dict(obj)
            % Calculate the value of G_dict based on G and keys
            value = array2table(obj.G);
            if ~isempty(obj.factorKeys)
                value.Properties.VariableNames = obj.factorKeys;
            end
            if ~isempty(obj.dimKeys)
                value.Properties.RowNames = obj.dimKeys;
            end
        end
        
        function value = get.A_dict(obj)
            % Calculate the value of A_dict based on A and keys
            if size(obj.A,2) ~= obj.nG
                A = zeros(0,obj.nG);
            else 
                A = obj.A;
            end
            value = array2table(A);
            if ~isempty(obj.factorKeys)
                value.Properties.VariableNames = obj.factorKeys;
            end
            if ~isempty(obj.conKeys)
                value.Properties.RowNames = obj.conKeys;
            end
        end

        function value = get.b_dict(obj)
            % Calculate the value of b_dict based on b and keys
            value = array2table(obj.b);
            if ~isempty(obj.conKeys)
                value.Properties.RowNames = obj.conKeys;
            end
        end

        function value = get.Z(obj)
            % Calculate conZono version
            % value = conZono;
            % value.G = obj.G;
            % value.c = obj.c;
            % value.A = obj.A;
            % value.b = obj.b;
            value = conZono(obj.G,obj.c,obj.A,obj.b);
        end

        function value = get.keys(obj)
            % Get struct of keys
            value.factors = obj.factorKeys;
            value.dims = obj.dimKeys;
            value.cons = obj.conKeys;
        end

        %% Setter Functions
        function obj = set.c_dict(obj, value)
            obj.c = table2array(value);
            obj.dimKeys = value.Properties.RowNames;
        end

        function obj = set.G_dict(obj, value)
            obj.G = table2array(value);
            obj.dimKeys = value.Properties.RowNames;
            obj.factorKeys = value.Properties.VariableNames;
        end
        
        function obj = set.A_dict(obj, value)
            obj.A = table2array(value);
            obj.conKeys = value.Properties.RowNames;
            obj.factorKeys = value.Properties.VariableNames;
        end

        function obj = set.b_dict(obj, value)
            obj.b = table2array(value);
        end

        function obj = set.Z(obj, value)
            % takes a ConZono and sets ConZonoM acordingly
            obj.c = value.c;
            obj.G = value.G;
            obj.A = value.A;
            obj.b = value.b;
        end

        function obj = set.keys(obj,value)
            % Sets the appropriete keys for a struct of keys
            if ~isempty(value.factors); obj.factorKeys = value.factors; end
            if ~isempty(value.dims); obj.dimKeys = value.dims; end
            if ~isempty(value.cons); obj.conKeys = value.cons; end
        end

        function obj = set.factorKeys(obj, value)
            obj.factorKeys = obj.setIndividualKeys(value,obj.nG);
        end

        function obj = set.dimKeys(obj, value)
            obj.dimKeys = obj.setIndividualKeys(value,obj.n);
        end

        function obj = set.conKeys(obj, value)
            obj.conKeys = obj.setIndividualKeys(value,obj.nC);
        end

        %% Basic Set Operations
        out = minkowskiSum(obj1,obj2)
        out = linearMap(obj1,obj2)
        out = generalizedIntersection(obj1,obj2,R,label)
        out = cartProd(obj1,obj2)

        % Define standard operations
        function out = plus(obj1,obj2)
            out = minkowskiSum(obj1,obj2);
        end
        function out = mtimes(obj1,obj2)
            out = linearMap(obj1,obj2);
        end

        function out = vertcat(varargin)
            out = varargin{1};%.copy();
            for i = 2:nargin
                out = cartProd(out,varargin{i});
            end                
        end
        
        function out = horzcat(varargin)
            out = varargin{1}.copy();
            for i = 2:nargin
                out = addInfo(out,varargin{i});
            end
        end
       
        %% Indexing
        % Overload Indexing
        B = subsref(A,S);   
        A = subsasgn(A,S,B); % untested and not really needed

        % Fancy Indexing (useless)
        function out = splitDims(obj,dims)
            S.type = '()'; S.subs = {dims};
            out = subsref(obj,S);
        end
        function out = splitFactors(obj,factors)
            S.type = '()'; S.subs = {':',factors};
            out = subsref(obj,S);
        end
        function out = splitCons(obj,cons)
            S.type = '()'; S.subs = {':',':',cons};
            out = subsref(obj,S);
        end
    
        %% Aditional Functions
        % Key Manipulation
        out = setIndividualKeys(obj,value,n)
        [keys1unique,keysBoth,keys2unique] = getUniqueFactorKeys(obj1,obj2)

        % Labeled Intersection
        out = labeledIntersection(obj1,obj2,dims,label)
                 
        % Add Info
        out = addInfo(obj1,obj2,dims,label);
    end
    

end