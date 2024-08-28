% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Class:
%       memZono
%       Dimension-aware and memory-encoded Zonotope
%       Z = {c + G \xi | ||\xi||_inf <= 1, A \xi = b, 
%               \xi_\in \{0,1\} \forall_{i \notin vset_}}
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
%       vset_ - nG x 1 logical vector defining if continous or discrete
%   Outputs:
%       Z - memZono object
%   Notes:
%       This class is built upon the functions written for the individual 
%       zonoLAB classes but adds dimensional awareness and memory-encoding/ 
%       preservation within the set operations.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
classdef memZono %< abstractZono %& matlab.mixin.CustomDisplay

    %% Data
    properties (Hidden,Access=private) % Underlying data structure
        G_      % Generator matrix
        c_      % Center
        A_ = [] % Constraint matrix
        b_ = [] % Constraint vector
        vset_    % vset_ defining if generators are continous or discrete
    end

    properties (Dependent,Hidden,Access=private) %(hidden/inaccessble)
        Gc_      % Continuous generator matrix (n x nGc)
        Gb_      % Binary generator matrix (n x nGb)
        Ac_      % Continuous constraint matrix (nC x nGc)
        Ab_      % Binary constraint matrix (nC x nGb)
    end

    properties (Dependent) % These properties get automatically updated when used
        c       % Center (n x 1)
        b       % Constraint vector (nC x 1)
        G       % Generator matrix (n x nG)
        Gc      % Continuous generator matrix (n x nGc)
        Gb      % Binary generator matrix (n x nGb)
        A       % Constraint matrix (nC x nG)
        Ac      % Continuous constraint matrix (nC x nGc)
        Ab      % Binary constraint matrix (nC x nGb)
        vset
    end

    % Dimensions
    properties (Dependent) % These properties get automatically updated when used
        n       % Dimension
        nG      % Number of generators
        nGc     % Number of continuous generators
        nGb     % Number of binary generators
        nC      % Number of constraints
    end
    properties (Dependent,Hidden)
        sz      % [n nG nC]
    end

    % I/O zono
    properties (Dependent,Hidden,Access=private) %<== unable to ensure private... ensure that dimensions are maintained...
        Z_          % Export to a base-zonotope class
    end
    properties (Dependent,Hidden)
        baseClass   % Equivalent base-zonotope class
    end

    % Labeling
    properties (Hidden,Access=protected)
        keys_ = struct( ...
            'factors',[],...
            'dims',[],...
            'cons',[])
    end
    properties (Dependent)
        keys
    end
    properties (Dependent,Hidden)
        factorKeys
        dimKeys
        conKeys
    end

    %% Constructors
    methods
        function obj = memZono(varargin)
            switch nargin
                case 1
                    if isa(varargin{1},'memZono') % <--- must have labels
                        obj = varargin{1};
                    else
                        error('A memZono must be created with labels.')
                    end
                case 2
                    obj.Z_ = varargin{1};
                    obj.keys = varargin{2};
                case 3
                    obj.Z_ = varargin{1};
                    obj.dimKeys = varargin{2};
                    obj.factorKeys = varargin{3};
                    obj.conKeys = varargin{3};
                case 4
                    obj.Z_ = varargin{1};
                    obj.dimKeys = varargin{2};
                    obj.factorKeys = varargin{3};
                    obj.conKeys = varargin{4};
                case 6
                    obj.G_ = varargin{1};
                    obj.c_ = varargin{2};
                    obj.A_ = varargin{3};
                    obj.b_ = varargin{4};
                    obj.vset_ = logical(varargin{5});
                    obj.keys = varargin{6};
                otherwise
                    error('Constructor not specified')
            end
        end
    end
    %% Parameter Set/Read
    methods
        % Matrices
        % Get Matrices
        function out = get.G_(obj) 
            if isempty(obj.G_); obj.G_ = zeros(obj.n,0); end
            out = obj.G_; 
        end
        function out = get.c_(obj); out = obj.c_; end %<== no checks included
        function out = get.A_(obj)
            if isempty(obj.A_); obj.A_ = zeros(0,obj.nG); end
            out = obj.A_; 
        end
        function out = get.b_(obj); out = obj.b_; end %<== no checks included

        % hybZono Matrices
        function out = get.Gc_(obj); out = obj.G_(:,obj.vset_); end
        function out = get.Gb_(obj); out = obj.G_(:,~obj.vset_); end
        function out = get.Ac_(obj); out = obj.A_(:,obj.vset_); end
        function out = get.Ab_(obj); out = obj.A_(:,~obj.vset_); end

        % Set hybZono Matrices
        function obj = set.Gc_(obj,in); obj.G_(:,obj.vset_) = in; end
        function obj = set.Gb_(obj,in); obj.G_(:,~obj.vset_) = in; end
        function obj = set.Ac_(obj,in); obj.A_(:,obj.vset_) = in; end
        function obj = set.Ab_(obj,in); obj.A_(:,~obj.vset_) = in; end

        % Dimensions
        function n = get.n(obj); n = size(obj.c_,1); end
        function nG = get.nG(obj); nG = size(obj.G_,2); end
        function nC = get.nC(obj); nC = size(obj.A_,1); end
        function sz = get.sz(obj); sz = [obj.n,obj.nG,obj.nC]; end
        % hybZono c/b factors
        function nGc = get.nGc(obj); nGc = sum(obj.vset_); end
        function nGb = get.nGb(obj); nGb = sum(~obj.vset_); end
    end

    %% In/Out with base Zonotope
    methods
        % test if special
        function out = issym(obj)
            % tests if any are symbolic
            if any([isa(obj.G_,'sym'),isa(obj.c_,'sym'),...
                    isa(obj.A_,'sym'),isa(obj.b_,'sym')])
                out = true;
            else
                out = false;
            end
        end
        function out = isnumeric(obj)
            % tests if all are numeric
            if all([isnumeric(obj.G_), isnumeric(obj.c_), ...
                    isnumeric(obj.A_), isnumeric(obj.b_)])
                out = true;
            else
                out = false;
            end
        end

        % Get base abstractZono class
        function out = get.baseClass(obj)
            if all(obj.vset_)
                if isempty(obj.A_), out = 'zono';
                else, out = 'conZono'; 
                end
            else, out = 'hybZono';
            end
        end

        % Output appropriate base zonotope
        function Z = get.Z_(obj)
            % warning('Ensure that your output dimensions line up correctly')
            switch obj.baseClass
                case 'zono'
                    Z = zono(obj.G_,obj.c_);
                case 'conZono'
                    Z = conZono(obj.G_,obj.c_,obj.A_,obj.b_);
                case 'hybZono'
                    Z = hybZono(obj.Gc_,obj.Gb_,obj.c_,...
                        obj.Ac_,obj.Ab_,obj.b_);
            end
        end

        % Output Specific Dimensions
        function out = Z(obj,dims)
            out = projection(obj,dims).Z_;
        end

        % Setter function for underlying zonotope data
        function obj = set.Z_(obj,in)
            switch class(in)
                case {'double','sym','optim.problemdef.OptimizationVariable'}
                    obj.c_ = in;
                    obj.vset_ = true(1,0);
                case 'zono'
                    obj.G_ = in.G;
                    obj.c_ = in.c;
                    obj.vset_ = true(1,in.nG);
                case 'conZono'
                    obj.G_ = in.G;
                    obj.c_ = in.c;
                    obj.A_ = in.A;
                    obj.b_ = in.b;
                    obj.vset_ = true(1,in.nG);
                case 'hybZono'
                    obj.G_ = [in.Gc,in.Gb];
                    obj.c_ = in.c;
                    obj.A_ = [in.Ac,in.Ab];
                    obj.b_ = in.b;
                    obj.vset_ = [true(1,in.nGc),false(1,in.nGb)];    
                case 'memZono'
                    obj.G_ = in.G_;
                    obj.c_ = in.c_;
                    obj.A_ = in.A_;
                    obj.b_ = in.b_;
                    obj.vset_ = in.vset_;              
            end
        end
    end
    
    %% Labeling 
    methods
        % Key Getter Functions
        function out = get.keys_(obj)
            out = obj.keys_;
        end
        function out = get.keys(obj)
            % out = obj.keys_; 
            out = structfun(@(x)(x'), obj.keys_,UniformOutput=false); %<== here for visibility in MATLAB
        end
        function out = get.dimKeys(obj); out = obj.keys_.dims; end
        function out = get.factorKeys(obj); out = obj.keys_.factors; end
        function out = get.conKeys(obj); out = obj.keys_.cons; end
            
        % Key Setter Functions
        function obj = set.keys(obj,in)
            % TODO: Add keys checks???
            if isstruct(in)
                obj.keys_ = in;
            else 
                obj.dimKeys = in;
                obj.factorKeys = in;
                obj.conKeys = in;
            end
        end
        function obj = set.dimKeys(obj,in)
            obj.keys_.dims = obj.keysCheck(in,obj.n);
        end
        function obj = set.factorKeys(obj,in)
            obj.keys_.factors = obj.keysCheck(in,obj.nG); 
        end
        function obj = set.conKeys(obj,in)
            obj.keys_.cons = obj.keysCheck(in,obj.nC);
        end
    end

    %% Generate Keys
    methods (Static)
        function [out] = genKeys(prefix,nums)
            labeler = @(prefix,num)sprintf('%s_%i',prefix,num);
            out = arrayfun(@(num){labeler(prefix,num)},nums);
        end
    end
    
    methods
        % Create keys from a patern
        function keys_ = keysStartsWith(obj,pattern)
            for field = string({'dims','factors','cons'})
                keys_.(field) = {};
                for i = 1:length(obj.keys_.(field))
                    if startsWith(obj.keys_.(field){i},pattern)
                        keys_.(field){end+1} = obj.keys_.(field){i};
                    end
                end
            end
        end
        
        % Relabel Dims
        function out = relabelDims(obj,inDims,outDims)
            if ~iscell(inDims); inDims = obj.keysStartsWith(inDims).dimKeys; end
            if ~iscell(outDims); outDims = memZono.genKeys(outDims,1:numel(inDims)); end
            out = obj.projection(inDims);
            out.dimKeys = outDims;
        end
    end
    %% Internal Keys Operations
    methods (Static)
        function out = keysCheck(in,n)
            % keysCheck(in,n) - checks to ensure the keys(in) is structured
            % currently for a dimension of n
            if n == 0; out = []; return; end
            if ~iscell(in); in = cellstr(in); end
            if length(in) == n
                out = in; 
            elseif isscalar(in)
                out{n} = [];
                for i = 1:n
                    out{i} = sprintf('%s_%d',in{1},i);
                end
            elseif length(in) ~= n
                error('keys not assigned correctly/wrong size');
            elseif numel(unique(out))<numel(out)
                error('Duplicate keys')
            else
                error('keys broken somehow');
            end
        end

        % Get all unique keys
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

        % Finds key indices
        function [idxk1,idxks1,idxks2,idxk2] = getKeyIndices(in1,in2)
            [k1,ks,k2] = memZono.getUniqueKeys(in1,in2);

            % Find Indices of individual keys
            idxk1 = zeros(1,length(k1));
            idxk2 = zeros(1,length(k2));
            for k = 1:length(k1), idxk1(k) = find(strcmp(in1,k1{k})); end
            for k = 1:length(k2), idxk2(k) = find(strcmp(in2,k2{k})); end

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

        function [lbl,idx] = getAllKeyIndices(obj1,obj2)
            %% Keys
            % shared dims
            [lbl.d1,lbl.ds,lbl.d2] = memZono.getUniqueKeys(obj1.dimKeys,obj2.dimKeys);
            [idx.d1,idx.ds1,idx.ds2,idx.d2] = memZono.getKeyIndices(obj1.dimKeys,obj2.dimKeys);
            % shared factors
            [lbl.k1,lbl.ks,lbl.k2] = memZono.getUniqueKeys(obj1.factorKeys,obj2.factorKeys);
            [idx.k1,idx.ks1,idx.ks2,idx.k2] = memZono.getKeyIndices(obj1.factorKeys,obj2.factorKeys);
            % shared cons
            [lbl.c1,lbl.cs,lbl.c2] = memZono.getUniqueKeys(obj1.conKeys,obj2.conKeys);
            [idx.c1,idx.cs1,idx.cs2,idx.c2] = memZono.getKeyIndices(obj1.conKeys,obj2.conKeys);
        end
    end

    %% General Methods
    methods
        %% Set Operations
        obj = map(obj1,obj2,inDims,outDims); % Maping function
        obj = transform(obj1,obj2,M,inDims,outDims); % Affine Mapping w/ dims (outdated/to be replaced)
        obj = memoryIntersection(obj1,obj2,sharedDimLabels); % Intersection
        obj = memorySum(obj1,obj2); % Minkowski Sum
        obj = memoryCartProd(obj1,obj2); % cartisian product
        obj = cartProd(obj1,obj2,dims1,dims2,options); % external cartisian product
    end

    %% Indexing  ----------------------------
    methods
        B = subsref(A,S);
        % A = subsasgn(A,S,B); %<---- not completed        

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
    end

    %% Ploting ----------------------------
    methods
        plot(obj,varargin);
    end

    %% Overloading ----------------------------
    methods
        % +, plus() - override for memorySum()
        function out = plus(in1,in2)
            out = in1.memorySum(in2);
        end
        % *, mtimes() - override for transform (map?) ... mainly just scalers
        function out = mtimes(in1,in2)
            if isa(in2,'memZono')
                out = in2.transform([],in1);%map(in2,in1); %<== flip the syntax order
            else
                error('mtimes only overloaded for one direction');
            end
        end
        % &, and() - overide memoryIntersection w/ checks
        function out = and(obj1,obj2,sharedDimLabels)
            if nargin == 2
                error('Must Supply labels for shared dimensions')
            end
            out = memoryIntersection(obj1,obj2,sharedDimLabels);
        end

        % Overide for memoryUnion w/ checks ... NOT IMPLIMENTED
        function obj = or(obj1,obj2)
            error('Union not yet implimented')
        end

        % vertcat (Extended cartProd)
        function obj = vertcat(varargin)
            warning('vertcat is not efficient yet')
            obj = varargin{1};
            for i = 2:nargin %<========= not efficient
                obj = cartProd(obj,varargin{i});
            end
        end

        % horzcat not yet decided (union? sum?) ... NOT IMPLIMENTED
        function obj = horzcat(varargin)
            if nargin == 1, obj = varargin{1};
            else, error('horzcat not defined');
            end
        end
    end

    methods (Static)
        % sum (extended plus)
        function obj = sum(varargin)
            obj = varargin{1};
            for i = 2:nargin %<========= not efficient
                obj = plus(obj,varargin{i});
            end
        end
        % all (extended and)
        function obj = all(varargin)
            obj = varargin{1};
            for i = 2:nargin %<========= not efficient
                obj = memoryIntersection(obj,varargin{i},sprintf('_all_%d',i));
            end
        end
        % any (extended or) ... NOT IMPLIMENTED
        function obj = any(varargin)
            obj = varargin{1};
            for i = 2:nargin
                obj = or(obj,varargin{2});
            end
        end
    end

    %% Input/Output and Display
    methods
        % zono data as tables
        function out = get.G(obj)
            out = array2table(obj.G_, RowNames=obj.dimKeys, VariableNames=obj.factorKeys); 
        end
        function out = get.c(obj)
            out = array2table(obj.c_, RowNames=obj.dimKeys, VariableNames={'c'}); 
        end
        function out = get.A(obj) 
            if obj.nC == 0, out = []; return; end
            out = array2table(obj.A_,RowNames=obj.conKeys,VariableNames=obj.factorKeys); 
        end
        function out = get.b(obj)
            if obj.nC == 0, out = []; return; end
            out = array2table(obj.b_,RowNames=obj.conKeys, VariableNames={'b'}); 
        end
        function out = get.vset(obj)
            out = array2table(reshape(obj.vset_,[],1),RowNames=obj.factorKeys,VariableNames={'vset_'});
        end

        %% Hybzono versions as tables
        function out = get.Gc(obj)
            out = obj.G(:,obj.keys.factors(obj.vset_));
        end
        function out = get.Gb(obj)
            if obj.nGb > 0, out = obj.G(:,obj.keys.factors(~obj.vset_));
            else, out = [];
            end
        end
        function out = get.Ac(obj)
            if obj.nC > 0, out =  obj.A(:,obj.keys.factors(obj.vset_));
            else, out = [];
            end
        end
        function out = get.Ab(obj) 
            if obj.nGb > 0 && obj.nC > 0
                out = obj.A(:,obj.keys.factors(~obj.vset_));
            else, out = [];
            end
        end
    end



    %%%%%%%%%%%%%%%%%%%
    % Questionable Usefulness/Not for Release
    %%%%%%%%%%%%%%%%%%%
    methods (Static)
        % Construction based on name of base objects
        %   (a helper method to combine multiple memZono/abstractZono objects based on name/size)
        function obj = varNameConstructor(varargin)
            for i = 1:nargin
                varName = inputname(i);
                obj_{i} = memZono(varargin{i},varName);
            end
            obj = vertcat(obj_{:});
        end
    end

    methods

        %% Affine (technically not needed but useful version instead of map version directly)
        function out = affine(in,M,b,inDims,outDims)
            out = in(inDims).transform(b,M,inDims,outDims);
            % out = in.map(M,inDims,outDims) + memZono(b,outDims);
        end

        % % Relabel all keys by adding a suffix
        % function out = relabel(obj,s,fields)
        %     arguments
        %         obj
        %         s
        %         fields = {'dims','factors','cons'}
        %     end
        %     keys_ = obj.keys_;
        %     for field = string(fields)
        %         for i = 1:numel(obj.keys_.(field))
        %             keys_.(field){i} = append(obj.keys_.(field){i},s);
        %         end
        %     end
        %     out = memZono(obj,keys_);
        % end

        %%% These are just direct versions of abstractZono methods... (not for initial release)
        % dimAwareFun
        varargout = dimAwareFun(fun,obj,dimIn,dimOut,lbl,options);
        % boundingBox
        function out = boundingBox(in); out = dimAwareFun(@boundingBox,in); end
        % convexHull
        function out = convexHull(in1,in2); out = dimAwareFun(@convexHull,in1,in2); end

    end
            


        % make these protected/private? --------------------------------------
    methods 
        % Output specific properties
        function varargout = varOut(obj,varargin)
            for i = 1:numel(varargin)
                varargout{i} = obj.(varargin{i});
            end
        end
        % output all data
        function varargout = exportAllData(obj)
            fields = {'G_','c_','A_','b_','vset_','keys_'};
            for i = 1:numel(fields); varargout{i} = obj.(fields{i}); end
        end

    end

















end

