classdef conTest < abstractTest

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

         function obj = conTest(G,c,A,b)
             obj.G = G;
             obj.c = c;
             obj.A = A;
             obj.b = b;
         end

        % get/set functions
        function nC = get.nC(obj); nC = size(obj.A,1); end

     end
end
