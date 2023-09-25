function out = generalizedIntersection(obj1,obj2,R,label)
% Generalized Intersection
arguments
    obj1    conZonoM
    obj2    conZonoM
    R       double
    label   char        = ''
end

% error checking
if (obj1.n ~= size(R,2)) || (obj2.n ~= size(R,1))
    error('Inconsistent dimensions!')
end

% constraint labels
if strcmp(label,'')
    conKeysNew = {};
else
    conKeysNew = constructKeys(label,(1:size(R,1)));
end

% note: need to fix the mess that is the A_dict assignment
[k1,kshared,k2] = getUniqueFactorKeys(obj1,obj2);
out = conZonoM;

out.c = obj1.c; % not dict version
out.G_dict = [obj1.G_dict(:,[k1,kshared]), zeroTable(obj1.n,k2)];
out.A_dict = [
    obj1.A_dict(:,[k1,kshared]), zeroTable(obj1.nC,k2); % obj1 constraints
    zeroTable(obj2.nC,k1), obj2.A_dict(:,[kshared,k2]); % obj2 constraints
    array2table(R*table2array(obj1.G_dict(:,k1)), ...
        "VariableNames",k1,"RowNames",conKeysNew), ...
        array2table( ...
            R*table2array(obj1.G_dict(:,kshared)) ...
            - table2array(obj2.G_dict(:,kshared)), ...
                "VariableNames", kshared,"RowNames",conKeysNew),...
        array2table(-1*table2array(obj2.G_dict(:,k2)), ...
            "VariableNames",k2,"RowNames",conKeysNew),
                % -1.*obj2.G_dict(:,k2) % intersetion constraint
];
out.b = [
    obj1.b; 
    obj2.b; 
    obj2.c-R*obj1.c
]; % not dict version


% all the labels are infered when seting the _dict versions
% out.factorKeys = [k1,kshared,k2];
% out.dimKeys = obj1.dimKeys;
% out.conKeys = [obj1.conKeys, obj2.conKeys, ...
%     constructKeys(label,1:size(R,1))]; 

% out.A_dict = [
%     obj1.A_dict(:,[k1,kshared]), zeroTable(obj1.nC,k2);
%     zeroTable(obj2.nC,k1), obj2.A_dict(:,[kshared,k2]);
%     array2table(R*table2array(obj1.G_dict(:,k1)),"VariableNames",k1), ...
%         array2table(R*table2array(obj1.G_dict(:,kshared)), ...
%             "VariableNames", kshared,"RowNames",conKeysNew) ...
%             + -1.*obj2.G_dict(:,kshared), ...
%         -1.*obj2.G_dict(:,k2)
% ];


end