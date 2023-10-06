classdef baseTest

    properties
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
    end

     methods

         % Constructors
         function obj = baseTest(G,c)
             obj.G = G;
             obj.c = c;
         end

     end
end
