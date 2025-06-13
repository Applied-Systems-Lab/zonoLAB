% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Generates a p x q matrix given an m x n kernel for convolutional type
%   operations
%           * implmentation only accounts for two dimensions
%   Method:
%       Returns a single matrix of p x q
%   Syntax:
%       [outMat]   = toplizify(weights,input)
%   Inputs:
%       weights - m x n kernel 
%       input -   r x s matrix
%   Outputs:
%       outMat - p x q matrix
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [outMat] = toplizify(weights,input)
    outMat = [];
    iH = size(input,1);
    iL = size(input,2);
    for w = 1: size(weights,4)
        tempw = weights(:,:,:,w);
        wH = size(weights,1);
        wL = size(weights,2);
        % 5 is relative to size of kernel
        for i = 1: wH
            % in actual implementation values should be fixed
                % the zeros padded (31) corresponds to input layer - size of kernel
            temp = triu(toeplitz([tempw(i,:) zeros(1,iL-wL)]));
                % 36 is the input layer
                % 32 is output
            temp = kron(diag(ones(iH,1),i-1), double(temp(1:iL-wL+1,1:iH)));
            out(:,:,i) = temp(1:(iL-wL+1)^2,1:iH^2);
        end
        outMat = [outMat; sum(out,3)];
    end
end

