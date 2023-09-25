function out = setIndividualKeys(~,value,n)
    % Update keys and ensure that it if a single value is passed
    % each generator is numbered distinctly
    if isempty(value)
        out = value;
    elseif isa(value,'cell')
        if length(value) == 1 && n ~= 1
            out = constructKeys(value,1:n);
        else
            out = value;
        end
    elseif isa(value,'string')
        out = constructKeys(value,1:n);
    else
        error('issues w/ constructor label')
    end
    
    out = reshape(out,1,[]);

end