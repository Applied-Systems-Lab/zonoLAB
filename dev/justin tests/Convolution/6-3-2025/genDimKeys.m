% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Generates Dimension Keys for a memZono object
%           * implmentation only accounts for two dimensions
%   Method:
%       Returns cell array of strings of form: [name]_[idx1]_[idx2]
%   Syntax:
%       [outKeys]   = genDimKeys(name, idx1, idx2)
%   Inputs:
%       name - title of these subset of keys
%       idx1 - 1 x n array counting along one dimension
%       idx2 - 1 x m array counting along another
%   Outputs:
%       outKeys - cell array of strings 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [outKeys] = genDimKeys(name, idx1, idx2)
    convKey = @(w,i,j) sprintf('%s_%d_%d',w,i,j);
    convKeys = @(w,row,col) arrayfun(@(i,j) convKey(w,i,j), row, col , UniformOutput=false);
    outKeys = convKeys(name,idx1,idx2)'; 
end

 