function keys = constructKeys(label, nums)
    % constructKeys: Constructs cell array of keys by concatenating a label
    % with each element of the nums array.
    %
    % Inputs:
    %   - label: The label string to be concatenated with the elements of nums.
    %   - nums: An array of numbers.
    %
    % Output:
    %   - keys: Cell array of keys constructed by concatenating the label
    %     with each element of nums.
    
    % keys = cellfun(@(x) [label, '_', num2str(x)], ...
    %     num2cell(nums), 'UniformOutput', false);
    if iscell(label); label = label{:}; end
    keys = {};
    for i = 1:length(nums)
        keys{i} = [label, '_', num2str(nums(i))];
    end
end
