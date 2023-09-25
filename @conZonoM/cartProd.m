function out = cartProd(obj1,obj2)
% Computes the cartesian product of two constrained zonotopes
% Inputs:
%   - obj1: conZonoM
%   - obj2: conZonoM
% Outputs:
%   - out: cartesian product (obj1 X obj2) as a conZonoM acounting for
%   common factors
%   
    arguments
       obj1
       obj2 conZonoM
    end

    % Acounting for factors
    [k1,kshared,k2] = getUniqueFactorKeys(obj1,obj2);
    out = conZonoM;
    out.factorKeys = [k1,kshared,k2];
    out.dimKeys = [obj1.dimKeys,obj2.dimKeys];
    out.conKeys = [obj1.conKeys,obj2.conKeys];

    % Performing set operations
    out.c = [
        obj1.c;
        obj2.c
    ];
    out.G_dict = [
        obj1.G_dict(:,k1), obj1.G_dict(:,kshared), zeroTable(obj1.n,k2);
        zeroTable(obj2.n,k1), obj2.G_dict(:,kshared), obj2.G_dict(:,k2)
    ];
    out.A_dict = [
        obj1.A_dict(:,k1), obj1.A_dict(:,kshared), zeroTable(obj1.nC,k2);
        zeroTable(obj2.nC,k1), obj2.A_dict(:,kshared), obj2.A_dict(:,k2)
    ];
    out.b = [
        obj1.b;
        obj2.b
    ];

end


