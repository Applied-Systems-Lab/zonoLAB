classdef conTest < baseTest

    properties
        A
        b
    end

    properties (Dependent)
        nC
    end


     methods

         function obj = conTest(G,c,A,b)
             obj@baseTest(G,c);
             obj.A = A;
             obj.b = b;
         end

        % get/set functions
        function nC = get.nC(obj); nC = size(obj.A,1); end

     end
end
