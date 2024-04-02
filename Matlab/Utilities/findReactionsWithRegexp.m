function result = findReactionsWithRegexp(model,pattern)       
    % Use regexp to find elements matching the pattern
    matches = regexp(model.rxns, pattern);
    
    % Convert matches cell array to logical array indicating matches
    matchesLogical = ~cellfun(@isempty, matches);
    
    % Use logical indexing to extract matched elements
    result = model.rxns(matchesLogical);
end