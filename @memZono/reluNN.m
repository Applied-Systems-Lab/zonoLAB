% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Remix for to return a memZono
%   Method:
%       Returns the memZono representation of a ReLU neural network (NN).
%   Syntax:
%       [NN]   = reluNN(X,Ws,bs,a)
%       [NN,Y] = reluNN(X,Ws,bs,a)
%   Inputs:
%       X  - zonotopic set in R^n (memeZono object)
%       Ws - 1 x nL cell array of weight matrices for the nL layers
%       bs - 1 x nL cell array of bias vectors for the nL layers
%       a  - positive scalar bounding the domain of each activation
%            function to the interval [-a,a]
%   Outputs:
%       NN - memZono in R^m representing NN function
%       Y  - memZono in R^m representing output set of NN function
%            corresponding to input set X
%   Notes:
%       Assumes all activation functions are ReLU. Assumes the input to
%       each activation function is within the interval [-a,a].
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [NN,Y] = reluNN(X,Ws,bs,a)
    arguments
        X memZono
        Ws (1,:) cell
        bs (1,:) cell
        a = 1000;
    end
    

    c = a/2*[0;1];
    Gc = a/2*[1,1,0,0;1,0,0,0];
    Gb = [0;0];
    Ac = a/2*[1,1,-1,-1;1,0,-1,0];
    Ab = a/2*[-2;-1];
    b = a/2*[0;-1];

    % Hybrid Zonotope of a single ReLU
    relu = hybZono(Gc,Gb,c,Ac,Ab,b);

    % Neural network memZono *assumes the pass in is a memZono and it is oriented corrently
    x0s = memZono.genKeys('x_L0_u',1:X.n);
    NN = relabelDims(X,X.dimKeys,x0s);

    prev_xs = x0s;
    for i = 1:(length(bs)-1)
        [n2,n1] = size(Ws{i});

        % v's are inputs to ReLU, x's are outputs of ReLU
        layer = memZono.empty();
        for j = 1:n2
            relu_ij = memZono(relu,{sprintf('v_L%d_u%d_',i,j),sprintf('x_L%d_u%d_',i,j)},sprintf('phi_L%d_u%d_',i,j));
            layer = cartProd(layer,relu_ij); % no new constraints will be added, so not providing labels
        end

        vs = layer.keysStartsWith('v').dims; % inputs to a ReLU layer
        xs = layer.keysStartsWith('x').dims; % outputs of a ReLU layer

        prev_layer = NN.map(Ws{i},prev_xs,vs)+memZono(bs{i},vs);
        layer = layer.and(prev_layer,sprintf('merge_L%i',i));
        NN = cartProd(NN,layer);
        prev_xs = xs;
    end

    % make connections from final hidden layer to output
    ys = memZono.genKeys('y',1:length(bs{end}));
    prev_layer = map(NN.projection(xs),Ws{end},xs,ys) + memZono(bs{end},ys);
    NN = cartProd(NN,prev_layer);

    NN = NN.projection([x0s,ys]); % input-output map
    Y = NN.projection(ys);  % output
end