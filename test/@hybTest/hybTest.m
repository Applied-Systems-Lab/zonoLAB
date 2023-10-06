classdef hybTest < conTest

    properties (Hidden)
        vset
    end

    properties (Dependent)
        Gc
        Gb
        Ac
        Ab
    end

    properties (Dependent)
        nGc
        nGb
    end


     methods

         function obj = hybTest(Gc,Gb,c,Ac,Ab,b)
             obj@conTest([Gc Gb],c,[Ac Ab],b);
             obj.vset = [repmat('c',size(Gc,2)),repmat('d',size(Gb,2))];
         end

        % get/set functions
        function Gc = get.Gc(obj)
            Gc = obj.G(:,obj.vset=='c');
        end
        function Gb = get.Gb(obj)
            Gb = obj.G(:,obj.vset=='d');
        end
        function Ac = get.Ac(obj)
            Ac = obj.A(:,obj.vset=='c');
        end
        function Ab = get.Ab(obj)
            Ab = obj.A(:,obj.vset=='d');
        end

        function nGc = get.nGc(obj); nGc = sum(obj.vset=='c'); end
        function nGb = get.nGb(obj); nGb = sum(obj.vset=='d'); end

     end
end
