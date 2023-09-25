function plotConZonoArray(Z)
    cmap = colormap(jet(length(Z)));
    % figure
    hold on
    for k = 1:length(Z)
        plot(Z{k}.Z,cmap(k,:), 0.5)
    end
    legend
    title(inputname(1))
end