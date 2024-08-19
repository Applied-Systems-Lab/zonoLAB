%% Plot Function for memZono
% - This calls the appropriate plot method for the 
%   base class with a specification of specific dimensions
function plot(obj,varargin)

    if obj.n > 3; error('Too many dimensions to plot'); end


    % Call plot command.
    Z_ = obj.Z_;
    plot(Z_,varargin{:});

    % Label Axis with dims
    xlabel(obj.dimKeys{1});
    ylabel(obj.dimKeys{2});
    if obj.n == 3; zlabel(obje.dimKeys{3}); end

end



% %% Old version
% function plot(obj,dims,varargin)
%     arguments
%         obj memZono
%         dims = [];
%     end
%     arguments (Repeating)
%         varargin
%     end
    
%     if isempty(dims), dims = obj.dimKeys;
%     elseif isnumeric(dims), dims = obj.dimKeys(dims);
%     end

%     Z_ = obj.Z(dims);
%     if Z_.n > 3, error('specify dims... too many to plot'); end
%     plot(Z_,varargin{:});
% end