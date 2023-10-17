classdef dataTest

    properties
        G
        c
        A
        b
    end

    properties (Dependent)
        n
        nG
        nC
    end

    methods
        function n = get.n(obj); n = size(obj.G,1); end
        function nG = get.nG(obj); nG = size(obj.G,2); end
        function nC = get.nC(obj); nC = size(obj.A,1); end
    end
end