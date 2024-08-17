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
    properties (Dependent,Hidden,Access=private) %<== ensure that dimensions are maintained...
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
        function out = get.G_(obj) 
            if isempty(obj.G_); obj.G_ = zeros(obj.n,0); end
            out = obj.G_; 
        end
        % function out = get.c_(obj); out = obj.c_; end
        function out = get.A_(obj)
            if isempty(obj.A_); obj.A_ = zeros(0,obj.nG); end
            out = obj.A_; 
        end
        % function out = get.b_(obj); out = obj.b_; end

        % % Set Matrices
        % function obj = set.G_(obj,in); obj.G_ = in; end
        % function obj = set.c_(obj,in); obj.c_ = in; end
        % function obj = set.A_(obj,in); obj.A_ = in; end
        % function obj = set.b_(obj,in); obj.b_ = in; end

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
        % hybZono dims
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

        function out = get.baseClass(obj)
            if all(obj.vset_)
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
        function out = get.keys(obj); 
            % out = obj.keys_; 
            out = structfun(@(x)(x'), obj.keys_,UniformOutput=false);
        end
        function out = get.dimKeys(obj); out = obj.keys_.dims; end
        function out = get.factorKeys(obj); out = obj.keys_.factors; end
        function out = get.conKeys(obj); out = obj.keys_.cons; end

        % Key Setter Functions
        function obj = set.keys(obj,in) %<-- add better checks?
            if isstruct(in)
                obj.keys_ = in;
            else 
                % error("don't set this way")
                obj.dimKeys = in;
                obj.factorKeys = in;
                obj.conKeys = in;
            end
        end
        function obj = set.dimKeys(obj,in)
            try obj.keys_.dims = obj.keysCheck(in,obj.n);
            catch; error('dim key set issue'); 
            end
        end
        function obj = set.factorKeys(obj,in)
            try obj.keys_.factors = obj.keysCheck(in,obj.nG); 
            catch; error('factor key set issue');
            end
        end
        function obj = set.conKeys(obj,in)
            try obj.keys_.cons = obj.keysCheck(in,obj.nC);
            catch; error('con key set issue');
            end
        end

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
            % out.dims = {};
            % for i=1:length(obj.dimKeys)
            %     if startsWith(obj.dimKeys{i},pattern)
            %         out.dims = [out.dims,obj.dimKeys{i}];
            %     end
            % end
            % out.factors = {};
            % for i=1:length(obj.factorKeys)
            %     if startsWith(obj.factorKeys{i},pattern)
            %         out.factors = [out.factors,obj.factorKeys{i}];
            %     end
            % end
            % out.cons = {};
            % for i=1:length(obj.conKeys)
            %     if startsWith(obj.conKeys{i},pattern)
            %         out.cons = [out.cons,obj.conKeys{i}];
            %     end
            % end
        end

        % Relabel all keys by adding a suffix
        function out = relabel(obj,s,fields)
            arguments
                obj
                s
                fields = {'dims','factors','cons'}
            end
            keys_ = obj.keys_;
            for field = string(fields)
                for i = 1:numel(obj.keys_.(field))
                    keys_.(field){i} = append(obj.keys_.(field){i},s);
                end
            end
            out = memZono(obj,keys_);
            % obj.dimKeys = memZono.genKeys(obj.dimKeys,s);
            % out = obj;
            % for i = 1: length(obj.dimKeys)
            %     obj.dimKeys{i} = append(obj.dimKeys{i}, s);
            % end
            % obj.factorKeys = memZono.genKeys(obj.factorKeys,s);
            % for i = 1: length(obj.factorKeys)
            %     obj.factorKeys{i} = append(obj.factorKeys{i}, s);
            % end
            % obj.conKeys = memZono.genKeys(obj.factorKeys,s);
            % for i = 1: length(obj.conKeys)
            %     obj.conKeys{i} = append(obj.conKeys{i}, s);
            % end
            % out = obj;
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
            elseif isscalar(in)
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
            labeler = @(prefix,num)sprintf('%s_%i',prefix,num);
            out = arrayfun(@(num){labeler(prefix,num)},nums);
            % if ~iscell(prefix)
            %     out = arrayfun(@(num){labeler(prefix,num)},nums);
            % else
            %     out = {};
            %     for i = 1:numel(prefix)
            %         out = [out, arrayfun(@(num){labeler(prefix{i},num)},nums)];
            %     end
            % end
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
            if ~iscell(inDims)
                if size(M,2) == 1; inDims = {inDims}; 
                else; inDims = memZono.genKeys(inDims,1:size(M,2)); end
            end
            if ~iscell(outDims)
                if size(M,1) == 1; outDims = {outDims}; 
                else;  outDims = memZono.genKeys(outDims,1:size(M,1)); end
            end
            out = in.transform([],M,inDims,outDims,retainExtraDims=false);%.projection(outDims);
        end

        % Copy constructor (allows relabeling dimension)
        function out = copy(obj,inDims,outDims)
            if ~iscell(inDims); inDims = obj.keysStartsWith(inDims).dimKeys; end
            if ~iscell(outDims); outDims = memZono.genKeys(outDims,1:numel(outDims)); end
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
        % plot(obj,dims,varargin);

        %% Overloading ----------------------------
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

        function obj = cartProd(obj1,obj2,dims1,dims2,s1,s2,options)
            arguments
                obj1 memZono
                obj2 memZono
                dims1 = [];
                dims2 = [];
                s1 = '_s1';
                s2 = '_s2';
                options.sharedMethod = 'reject';
            end
            if isempty(dims1); dims1 = obj1.dimKeys; end
            if isempty(dims2); dims2 = obj2.dimKeys; end
            [k1,ks,k2] = memZono.getUniqueKeys(dims1,dims2);
            switch options.sharedMethod
                case 'reject'
                    if ~isempty(ks)%intersect(dims1,dims2))
                        error('standard cartProd only works if no dims are in common')
                    end
                    obj = merge(obj1,obj2);
                case 'rename'
                    if ~isempty(ks)
                        Z1 = obj1.projection(k1); 
                        Z2 = obj2.projection(k2);
                        Z1s = obj1.projection(ks).relabel(s1,'dims');
                        Z2s = obj2.projection(ks).relabel(s2,'dims');
                        % perform cartProd()
                        obj = merge(Z1s,Z2s);
                        if Z1.n>0; obj.merge(Z1); end
                        if Z2.n>0; obj.merge(Z2); end
                    else
                        obj = merge(Z1,Z2);
                    end
                    % obj = merge(merge(Z1,Z1s),merge(Z2s,Z2));
                % case 'merge'
                %     obj = merge(obj1,obj2);
                % case 'combine'
                %     obj = combine(obj1,obj2,lbl_new);
            end
            warning('cartProd not implimented directly (uses merge w/ a check)... should reimpliment as vertcat and then overload cartProd() as vertcat');
        end

        % Extended intersection
        function obj = vertcat(varargin)
            warning('vertcat is not efficient yet')
            obj = varargin{1};
            for i = 2:nargin %<========= not efficient
                % obj = merge(obj,varargin{i});
                obj = cartProd(obj,varargin{i});
            end
        end

        % Extended minkowsi sum
        function obj = horzcat(varargin)
            error('horzcat not defined')
            % obj = varargin{1};
            % for i = 2:nargin %<========= not efficient
            %     obj = combine(obj,varargin{i});
            % end
        end

        %% Indexing
        B = subsref(A,S);
        % A = subsasgn(A,S,B); %<---- not completed
        

        % Projection is defined for internal use - subsref (indexing) is simpilar syntax
        function out = projection(obj,dims)
            if ~iscell(dims) % if not already in cell form
                dims = obj.keysStartsWith(dims).dimKeys;
            end
            [~,idx] = ismember(dims,obj.dimKeys);
            keys_ = obj.keys_; keys_.dims = dims;
            out = memZono(obj.G_(idx,:),obj.c_(idx,:),obj.A_,obj.b_,obj.vset_,keys_);
        end

    end

    %% Input/Output and Display
    % properties (Dependent)
    %     factorKeys_
    %     dimKeys_
    %     conKeys_
    % end

    methods
        function out = get.G(obj)
            out = array2table(obj.G_, RowNames=obj.dimKeys, VariableNames=obj.factorKeys); 
        end
        function out = get.c(obj)
            out = array2table(obj.c_, RowNames=obj.dimKeys, VariableNames={'c'}); 
        end
        function out = get.A(obj) 
            out = array2table(obj.A_,RowNames=obj.conKeys,VariableNames=obj.factorKeys); 
        end
        function out = get.b(obj)
            out = array2table(obj.b_,RowNames=obj.conKeys, VariableNames={'b'}); 
        end
        function out = get.vset(obj)
            out = array2table(reshape(obj.vset_,[],1),RowNames=obj.factorKeys,VariableNames={'vset_'});
        end

        function out = get.Gc(obj)
            out = array2table(obj.Gc_, RowNames=obj.dimKeys, VariableNames=obj.factorKeys(obj.vset_)); 
        end
        function out = get.Gb(obj)
            out = array2table(obj.Gb_, RowNames=obj.dimKeys, VariableNames=obj.factorKeys(~obj.vset_)); 
        end
        function out = get.Ac(obj) 
            out = array2table(obj.Ac_,RowNames=obj.conKeys,VariableNames=obj.factorKeys(obj.vset_)); 
        end
        function out = get.Ab(obj) 
            out = array2table(obj.Ab_,RowNames=obj.conKeys,VariableNames=obj.factorKeys(~obj.vset_)); 
        end

        % function out = get.factorKeys_(obj); out = reshape(obj.factorKeys,[],1); end
        % function out = get.dimKeys_(obj); out = reshape(obj.dimKeys,[],1); end
        % function out = get.conKeys_(obj); out = reshape(obj.conKeys,[],1); end

    end






end

