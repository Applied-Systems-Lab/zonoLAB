function out = labeledIntersection(obj1,obj2,dims,label)
% labeledIntersection - performs an intersection upon specific dims
arguments
    obj1    conZonoM
    obj2    conZonoM
    dims
    label   char = ''
end

if obj2.n > length(dims) %if correct size, obj2 labels are ignored
    subs.type = '()';
    subs.subs = dims;
    obj2 = subsref(obj2,subs); % grab only specified dims... (eliminates need for zero R rows)
end

% R is constructed to only operate on the dims selected from obj2
R = zeros(length(dims),obj1.n);
for i = 1:length(dims)
    j = strcmp(obj1.dimKeys,dims(i)); % numerical index of dim in obj1
    R(i,j) = 1;
end

out = generalizedIntersection(obj1,obj2,R,label);

end
