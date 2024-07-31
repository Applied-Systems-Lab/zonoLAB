% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Class:
%       memZono
%       Dimension-aware and memory-encoded Zonotope
%       Z = {c + G \xi | ||\xi||_inf <= 1, A \xi = b, 
%               \xi_\in \{0,1\} \forall_{i \notin vset}}
%       TODO: finalize definition @jruths - what else do we want for this
%   Syntax:
%       TODO: add
%   Inputs:
%       z - zono object in base zonotope form
%       keys - a struct specifying the underying keys or a string to generate keys from
%       G - n x nG matrix to define zonotope in R^n with nG generators
%       c - n x 1 vector to define center
%       A - nC x nG matrix to define nC equality constraints (A \xi = b)
%       b - nC x 1 vector to define nC equality constraints (A \xi = b)
%       vSet - nG x 1 logical vector defining if continous or discrete
%   Outputs:
%       Z - memZono object
%   Notes:
%       This class is built upon the functions written for the individual 
%       zonoLAB classes but adds dimensional awareness and memory-encoding/ 
%       preservation within the set operations.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
classdef memZono

    %% Data
    properties (Hidden) % Underlying data structure
        G_      % Generator matrix
        c_      % Center
        A_ = [] % Constraint matrix
        b_ = [] % Constraint vector
        vset    % vSet defining if generators are continous or discrete
    end

    properties (Dependent)
        c       % Center (n x 1)
        b       % Constraint vector (nC x 1)
        % G       % Generator matrix (n x nG)
        Gc      % Continuous generator matrix (n x nGc)
        Gb      % Binary generator matrix (n x nGb)
        % A       % Constraint matrix (nC x nG)
        Ac      % Continuous constraint matrix (nC x nGc)
        Ab      % Binary constraint matrix (nC x nGb)
    end

    properties (Dependent,Hidden) %(hidden do to hyb-zono precidence)
        G       % Generator matrix (n x nG)
        A       % Constraint matrix (nC x nG)
    end

    % Dimensions
    properties (Dependent) % These properties get automatically updated when used
        n       % Dimension
        nG      % Number of generators
        nGc     % Number of continuous generators
        nGb     % Number of binary generators
        nC      % Number of constraints
    end

    % I/O zono
    properties (Dependent, Hidden)
        Z           % Export to a base-zonotope class
        baseClass   % Equivalent base-zonotope class
    end

    % Labeling
    properties (Hidden)
        keys = struct( ...
            'factors',[],...
            'dims',[],...
            'cons',[])
    end
    properties (Dependent)
        factorKeys
        dimKeys
        conKeys
    end

    %% Constructors
    methods
        function obj = memZono(varargin)
            if nargin == 1
                if isa(varargin{1},'memZono') % <--- must have labels
                    obj = varargin{1};
                else
                    error('A memZono must be created with labels.')
                end
            elseif nargin == 2 % <--- base object and labels
                obj.Z = varargin{1};
                obj.keys = varargin{2};
            elseif nargin == 6 % <--- direct definition (primarily internal use)
                obj.G_ = varargin{1};
                obj.c_ = varargin{2};
                obj.A_ = varargin{3};
                obj.b_ = varargin{4};
                obj.vset = logical(varargin{5});
                obj.keys = varargin{6};
            else
                error('non-simple constructor not specified')
            end
        end
    end
    methods (Static)
        % Costruction based on name of base objects
        function obj = varNameConstructor(varargin)
            for i = 1:nargin
                varName = inputname(i);
                obj_{i} = memZono(varargin{i},varName);
            end
            obj = vertcat(obj_{:});
        end

    end

    %% Parameter Set/Read
    methods
        % Matrices
        % Get Matrices
        function out = get.G(obj) 
            if isempty(obj.G_); obj.G_ = zeros(obj.n,0); end
            out = obj.G_; 
        end
        function out = get.c(obj); out = obj.c_; end
        function out = get.A(obj)
            if isempty(obj.A_); obj.A_ = zeros(0,obj.nG); end
            out = obj.A_; 
        end
        function out = get.b(obj); out = obj.b_; end

        % Set Matrices
        function obj = set.G(obj,in); obj.G_ = in; end
        function obj = set.c(obj,in); obj.c_ = in; end
        function obj = set.A(obj,in); obj.A_ = in; end
        function obj = set.b(obj,in); obj.b_ = in; end

        % hybZono Matrices
        function out = get.Gc(obj); out = obj.G(:,obj.vset); end
        function out = get.Gb(obj); out = obj.G(:,~obj.vset); end
        function out = get.Ac(obj); out = obj.A(:,obj.vset); end
        function out = get.Ab(obj); out = obj.A(:,~obj.vset); end

        % Set hybZono Matrices
        function obj = set.Gc(obj,in); obj.G_(:,obj.vset) = in; end
        function obj = set.Gb(obj,in); obj.G_(:,~obj.vset) = in; end
        function obj = set.Ac(obj,in); obj.A_(:,obj.vset) = in; end
        function obj = set.Ab(obj,in); obj.A_(:,~obj.vset) = in; end

        % Dimensions
        function n = get.n(obj); n = size(obj.c,1); end
        function nG = get.nG(obj); nG = size(obj.G,2); end
        function nC = get.nC(obj); nC = size(obj.A,1); end
        % hybZono dims
        function nGc = get.nGc(obj); nGc = sum(obj.vset); end
        function nGb = get.nGb(obj); nGb = sum(~obj.vset); end


        % Additional Propterties
        function out = dimMin(obj,dims)
            out = projection(obj,dims).Z.lb;
        end
        function out = dimMax(obj,dims)
            out = projection(obj,dims).Z.ub;
        end
        function out = dimBounds(obj,dims)
            out = projection(obj,dims).Z.bounds;
        end
    end

    %% In/Out with base Zonotope
    methods
        % test if special
        function out = issym(obj)
            % tests if any are symbolic
            if any([isa(obj.G,'sym'),isa(obj.c,'sym'),...
                    isa(obj.A,'sym'),isa(obj.b,'sym')])
                out = true;
            else
                out = false;
            end
        end
        function out = isnumeric(obj)
            % tests if all are numeric
            if all([isnumeric(obj.G), isnumeric(obj.c), ...
                    isnumeric(obj.A), isnumeric(obj.b)])
                out = true;
            else
                out = false;
            end
        end

        function out = get.baseClass(obj)
            if all(obj.vset)
                if isempty(obj.A_)
                    out = 'zono';
                else
                    out = 'conZono';
                end
            else
                out = 'hybZono';
            end
        end

        % Output appropriate base zonotope
        function Z = get.Z(obj)
            % warning('Ensure that your output dimensions line up correctly')
            switch obj.baseClass
                case 'zono'
                    Z = zono(obj.G,obj.c);
                case 'conZono'
                    Z = conZono(obj.G,obj.c,obj.A,obj.b);
                case 'hybZono'
                    Z = hybZono(obj.Gc,obj.Gb,obj.c,...
                        obj.Ac,obj.Ab,obj.b);
            end
        end

        % Setter function for underlying zonotope data
        function obj = set.Z(obj,in)
            switch class(in)
                case {'double','sym','optim.problemdef.OptimizationVariable'}
                    obj.c_ = in;
                    obj.vset = true(1,0);
                case 'zono'
                    obj.G_ = in.G;
                    obj.c_ = in.c;
                    obj.vset = true(1,in.nG);
                case 'conZono'
                    obj.G_ = in.G;
                    obj.c_ = in.c;
                    obj.A_ = in.A;
                    obj.b_ = in.b;
                    obj.vset = true(1,in.nG);
                case 'hybZono'
                    obj.G_ = [in.Gc,in.Gb];
                    obj.c_ = in.c;
                    obj.A_ = [in.Ac,in.Ab];
                    obj.b_ = in.b;
                    obj.vset = [true(1,in.nGc),false(1,in.nGb)];    
                case 'memZono'
                    obj.G_ = in.G;
                    obj.c_ = in.c;
                    obj.A_ = in.A;
                    obj.b_ = in.b;
                    obj.vset = in.vset;              
            end
        end

    end
    
    %% Labeling
    methods
        % Key Getter Functions
        function out = get.keys(obj); out = obj.keys; end
        function out = get.factorKeys(obj); out = obj.keys.factors; end
        function out = get.dimKeys(obj); out = obj.keys.dims; end
        function out = get.conKeys(obj); out = obj.keys.cons; end

        % Key Setter Functions
        function obj = set.keys(obj,in)
            if isstruct(in) %<-- add better checks?
                obj.keys = in;
            else
                obj.factorKeys = in;
                obj.dimKeys = in;
                obj.conKeys = in;
            end
        end
        function obj = set.factorKeys(obj,in)
            try obj.keys.factors = obj.keysCheck(in,obj.nG); 
            catch
                warning('factor key set issue (toc used as key)');
                obj.keys.factors = obj.keysCheck(num2str(toc),obj.nG);
            end
        end
        function obj = set.dimKeys(obj,in)
            try obj.keys.dims = obj.keysCheck(in,obj.n);
            catch; error('dim key set issue'); 
            end
        end
        function obj = set.conKeys(obj,in)
            try obj.keys.cons = obj.keysCheck(in,obj.nC);
            catch
                warning('con key set issue (toc used as key)');
                obj.keys.cons = obj.keysCheck(num2str(toc),obj.nC);
            end
        end

        % Create keys from a patern
        % function out = cellStartsWith(in,pattern)
        %     idx = cellfun(@(lbl) startsWith(lbl,pattern),in);
        %     out = in(idx);
        %     % out = in();
        % end
        % function out = factorKeyStartsWith(obj,pattern)
        %     out = cellStartsWith(obj.factorKeys,pattern);
        % end
        % function out = dimKeyStartsWith(obj,pattern)
        %     out = cellStartsWith(obj.dimKeys,pattern);
        % end
        % function out = conKeysStartsWith(obj,pattern)
        %     out = cellStartsWith(obj.conKeys,pattern);
        % end
        function out = keysStartsWith(obj,pattern)
            % out.factorKeys = cellStartsWith(obj.factorKeys,pattern);
            % out.dimKey = cellStartsWith(obj.dimKeys,pattern);
            % out.conKeys = cellStartsWith(obj.conKeys,pattern);
            out.factorKeys = {};
            out.dimKeys = {};
            out.conKeys = {};
            for i=1:length(obj.factorKeys)
                if startsWith(obj.factorKeys{i},pattern)
                    out.factorKeys = [out.factorKeys,obj.factorKeys{i}];
                end
            end
            for i=1:length(obj.dimKeys)
                if startsWith(obj.dimKeys{i},pattern)
                    out.dimKeys = [out.dimKeys,obj.dimKeys{i}];
                end
            end
            for i=1:length(obj.conKeys)
                if startsWith(obj.conKeys{i},pattern)
                    out.conKeys = [out.conKeys,obj.conKeys{i}];
                end
            end
        end

        % Relabel all keys by adding a suffix
        function out = relabel(obj,s)
            for i = 1: length(obj.factorKeys)
                   obj.factorKeys{i} = append(obj.factorKeys{i}, s);
            end
            for i = 1: length(obj.dimKeys)
                   obj.dimKeys{i} = append(obj.dimKeys{i}, s);
            end
            for i = 1: length(obj.conKeys)
                   obj.conKeys{i} = append(obj.conKeys{i}, s);
            end
            out = obj;
        end
    end

    %% Keys Stuff
    methods (Static)
        function out = keysCheck(in,n)
            % keysCheck(in,n) - checks to ensure the keys(in) is structured
            % currently for a dimension of n
            if n == 0; out = []; return; end
            if ~iscell(in); in = cellstr(in); end
            if length(in) == n
                out = in; 
            elseif length(in) == 1
                out{n} = [];
                for i = 1:n
                    out{i} = sprintf('%s_%d',in{1},i);
                end
            elseif length(in) ~= n
                error('keys not assigned correctly/wrong size');
            else
                error('keys broken');
            end
        end

        function [k1,ks,k2] = getUniqueKeys(in1,in2)
            if isempty(in1) || isempty(in2)
                ks = {};
                k1 = in1;
                k2 = in2;
            else
                ks = intersect(in1,in2);
                k1 = setdiff(in1,ks);
                k2 = setdiff(in2,ks);
            end
        end

        function [idxk1,idxks1,idxks2,idxk2] = getKeyIndices(in1,in2)
            [k1,ks,k2] = memZono.getUniqueKeys(in1,in2);

            % Find Indices of individual keys
            idxk1 = zeros(1,length(k1));
            idxk2 = zeros(1,length(k2));
            for k = 1:length(k1)
                idxk1(k) = find(strcmp(in1,k1{k}));
            end
            for k = 1:length(k2)
                idxk2(k) = find(strcmp(in2,k2{k}));
            end


            % Find indices of shared keys
            if isempty(ks)
                idxks1 = []; idxks2 = [];
            else
                idxks1 = zeros(1,length(ks));
                idxks2 = zeros(1,length(ks));
                for k = 1:length(ks)
                    idxks1(k) = find(strcmp(in1,ks{k}));
                    idxks2(k) = find(strcmp(in2,ks{k}));
                end
            end
        end

        function [out] = genKeys(prefix,nums)
            labeler = @(letter,num)sprintf('%s_%i',letter,num);
            if ~iscell(prefix)
                out = arrayfun(@(num){labeler(prefix,num)},nums);
            else
                out = {};
                for i = 1:numel(prefix)
                    out = [out,arrayfun(@(num){labeler(prefix{i},num)},nums)];
                end
            end
        end

    end

    %% General Methods
    methods
        %% Set Operations
        obj = transform(obj1,obj2,M,inDims,outDims); % Affine Mapping w/ dims
        obj = merge(obj1,obj2,sharedDimLabels); % Intersection
        obj = combine(obj1,obj2); % Minkowski Sum

        % Additional Methods
        function out = linMap(in,M,inDims,outDims)
            % if isscalar(varargin)
            %     inDims = in.dimKeys;
            %     outDims = varargin{1}; 
            %     warning('lack of inDims specification can cause issues with dimension ordering'); 
            % else
            %     inDims = varargin{1};
            %     outDims = varargin{2};
            % end
            if ~iscell(inDims); inDims = memZono.genKeys(inDims,1:size(M,1)); end
            if ~iscell(outDims); outDims = memZono.genKeys(outDims,1:size(M,1)); end
            out = in.transform([],M,inDims,outDims,retainExtraDims=false);%.projection(outDims);
        end

        % Copy constructor (allows relabeling dimension)
        function out = copy(obj,inDims,outDims)
            if ~iscell(inDims); inDims = obj.keysStartsWith(inDims).dimKeys; end
            if ~iscell(outDims); outDims = memZono.genKeys(outDims,numel(outDims)); end
            out = obj.projection(inDims);
            out.dimKeys = outDims;
            % if nargin == 2
            %     warning('relabeling dimensions without specifying order (should be depreciated)')
            %     out.dimKeys = varargin{1};
            % elseif nargin == 3
            %     out = obj.projection(varargin{1});
            %     out.dimKeys = varargin{2};
            % end
        end
        
        %% Ploting
        plot(obj,dims,varargin);

        %% Overloading
        function out = plus(in1,in2)
            out = in1.combine(in2);
        end
        function out = mtimes(in1,in2)
            if isa(in2,'memZono')
                out = in2.transform([],in1); %<== flip the syntax order
            else
                error('mtimes only overloaded for direction')
            end
        end
        function out = and(obj1,obj2)
            out = merge(obj1,obj2);
        end

        function obj = or(obj1,obj2)
            error('Union not yet implimented')
            % obj = union(obj1,obj2);
        end
                
        % Extended intersection
        function obj = vertcat(varargin)
            obj = varargin{1};
            for i = 2:nargin %<========= not efficient
                obj = merge(obj,varargin{i});
            end
        end

        % Extended minkowsi sum
        function obj = horzcat(varargin)
            obj = varargin{1};
            for i = 2:nargin %<========= not efficient
                obj = combine(obj,varargin{i});
            end
        end

        %% Indexing
        B = subsref(A,S);
        % A = subsasgn(A,S,B); %<---- not completed
        

        % Projection is defined for internal use - subsref (indexing) is simpilar syntax
        function out = projection(obj,dims)
            if ~iscell(dims) % if not already in cell form
                if strcmp(dims,'all')
                    dims = obj.dimKeys; 
                else
                    dims = obj.keysStartsWith(dims).dimKeys;
                end
            end
            [~,idx] = ismember(dims,obj.dimKeys);
            keys_ = obj.keys;
            keys_.dims = dims;
            out = memZono(obj.G(idx,:),obj.c(idx,:),obj.A,obj.b,obj.vset,keys_);
        end

    end
end
