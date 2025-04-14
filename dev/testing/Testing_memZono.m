clear
Z = zono([eye(2),ones(2)],ones(2,1));
Z2 = zono([eye(4),2*eye(4)],zeros(4,1));
Z3 = union(cartProd(Z,Z),Z2);

temp = memZono(Z,'test');

dim_ = @(layer,i,j) genMatIndexKeys(sprintf("layer%d",layer), i,j)


% keys(1:2,4:5)

% temp2 = memZono()


M = sparse([0, 0, 2, 4; 1, 0, 0, 0]);
temp = memZono(Z2,'z2');
out = linMap(temp,M,'z2','z3')






% temp.dimKeys = 'test';
% temp2.dimKeys = {'test_2','test_3'};

% minSum(temp,temp2)
% temp+temp2

% temp3 = and(temp,temp2)

% plot(temp3)


function keys = genMatIndexKeys(key,I,J)
    keys = {};
    for i = I
        keys = [keys,memZono.genKeys(sprintf("%s_%d",key,i),J)];
    end
end