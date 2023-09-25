function out = linearMap(obj1,obj2)
arguments
    obj1 
    obj2 conZonoM
end

% Use ConZono Function
out = conZonoM;
out.Z = obj1 * obj2.Z; % uses the conZono class definition

% Setting Keys
out.keys = obj2.keys;


end