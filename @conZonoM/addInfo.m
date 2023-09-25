function out = addInfo(obj1,obj2,dims,conLabels)
% adds information to a zonotope... if 
arguments
    obj1    
    obj2
    dims    = {};
    conLabels   = '';
end

% input conditioning
if isempty(dims); dims = obj2.dimKeys; end
% if isempty(conLabels); conLabels = num2str(obj1.nC + obj2.nC + 1); end

% setup
S.type = '()';
S.subs = {dims};
obj2 = subsref(obj2,S);
[k1,kshared,k2] = getUniqueDimKeys(obj1,obj2);

S.subs = {[k1,kshared]};
out = subsref(obj1,S);

% Cartisian Product
if ~isempty(k2)
    S.subs = {k2};
    out = cartProd(out,subsref(obj2,S));
end

% Labeled Intersection
if ~isempty(kshared)
out = labeledIntersection(out,obj2,kshared,conLabels);

end
