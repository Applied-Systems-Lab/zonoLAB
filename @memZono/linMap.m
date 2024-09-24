function out = linMap(in,M,inDims,outDims)
    % Run linMap functionality
    if ~iscell(inDims)
        if size(M,2) == 1; inDims = {inDims}; 
        else; inDims = memZono.genKeys(inDims,1:size(M,2)); end
    end
    if ~iscell(outDims)
        if size(M,1) == 1; outDims = {outDims}; 
        else;  outDims = memZono.genKeys(outDims,1:size(M,1)); end
    end
    out = in.transform([],M,inDims,outDims,retainExtraDims=false);
end