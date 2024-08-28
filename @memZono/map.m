function out = map(obj1,obj2,inDims,outDims)
    % check to ensure only 
    memZonoBools = [isa(obj1,'memZono'),isa(obj2,'memZono')];
    if all(memZonoBools)
        error('Map currently not yet implimented for functional mapping');
    else
        if memZonoBools(1)
            in = obj1;
            M = obj2;
        elseif memZonoBools(2)
            in = obj2;
            M = obj1;
        else
            error('this should never be reached');
        end
        out = linMap(in,M,inDims,outDims);
    end

    % switch memZonoBools
    %     case [1 1]
    %         error('Map currently not yet implimented for functional mapping');
    %     case [1 0]
    %         in = obj1;
    %         M = obj2;
    %         out = linMap(in,M,inDims,outDims);
    %     case [0 1]
    %         in = obj2;
    %         M = obj1;
    %         out = linMap(in,M,inDims,outDims);
    %     case [0 0]
    %         error("this shouldn't be possible")
    % end
end

% Local linMap... since not a method currently
function out = linMap(in,M,inDims,outDims);
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

%     % Run linMap functionality
%     if ~iscell(inDims)
%         if size(obj2,2) == 1; inDims = {inDims}; 
%         else; inDims = memZono.genKeys(inDims,1:size(obj2,2)); end
%     end
%     if ~iscell(outDims)
%         if size(obj2,1) == 1; outDims = {outDims}; 
%         else;  outDims = memZono.genKeys(outDims,1:size(obj2,1)); end
%     end
%     out = obj1.transform([],obj2,inDims,outDims);


%     if all(memZonoBools)
%         error('Map currently not yet implimented for functional mapping');
%     else
        
%         % run linMap
%         in = obj1;
%         M = obj2;
%         if ~iscell(inDims)
%             if size(obj2,2) == 1; inDims = {inDims}; 
%             else; inDims = memZono.genKeys(inDims,1:size(obj2,2)); end
%         end
%         if ~iscell(outDims)
%             if size(obj2,1) == 1; outDims = {outDims}; 
%             else;  outDims = memZono.genKeys(outDims,1:size(obj2,1)); end
%         end
%         out = obj1.transform([],obj2,inDims,outDims);
%     end
% end