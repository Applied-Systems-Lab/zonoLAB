function t = zeroTable(n,keys)
    % creates a table of zeros for a given set of keys
    t = array2table(zeros(n,length(keys)),"VariableNames",keys);
end