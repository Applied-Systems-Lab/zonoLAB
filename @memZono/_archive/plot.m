%% Plot Function for memZono
% - This calls the appropriate plot method for the 
%   base class with a specification of specific dimensions
function plot(obj,varargin)

    if obj.n <= 0; error('Need dimensions to plot'); end
    if obj.n > 3; error('Too many dimensions to plot'); end

    % Call plot command.
    Z_ = obj.Z_;
    plot(Z_,varargin{:});

    % Label Axis with dims
    ax = gca;
    if isempty(ax.XLabel.String), xlabel(obj.dimKeys{1}); end
    if isempty(ax.YLabel.String), ylabel(obj.dimKeys{2}); end
    if obj.n == 3
        if isempty(ax.ZLabel.String), zlabel(obj.dimKeys{3}); end; 
    end

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