classdef abstractTest

    properties (Abstract)
        G
        c
    end

    properties (Dependent)
        n
        nG
    end

    methods
        % operations
        obj = plus(obj1,obj2)

        % get/set functions
        function n = get.n(obj); n = size(obj.G,1); end
        function nG = get.nG(obj); nG = size(obj.G,2); end
        function nC = get.nC(obj); nC = size(obj.A,1); end
    end

end