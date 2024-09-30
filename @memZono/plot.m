%% Plot Function for memZono
% - This calls the appropriate plot method for the 
%   base class with a specification of specific dimensions
function plot(obj,dims,varargin)
    % Project acording to dims    
    Z_ = obj.Z(dims);
    if Z_.n > 3, error('specify dims... too many to plot'); end
    % Call baseZono plot() method
    plot(Z_,varargin{:});
end