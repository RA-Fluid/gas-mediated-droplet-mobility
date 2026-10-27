function [rootsFound, rootStates] = RA_BGK_SolveBrackets(brackets, fixed, options, evaluator)
    % RA_BGK_SOLVEBRACKETS Solve all sign-changing thickness intervals.
    n = size(brackets, 1);
    rootsFound = NaN(n, 1);
    rootStates = repmat(RA_BGK_EmptyState(), n, 1);
    for k = 1:n
        rootsFound(k) = fzero(@residualOnly, brackets(k, :), ...
                              optimset('TolX', options.tolDelta, 'Display', 'off'));
        rootStates(k) = evaluator(rootsFound(k), fixed, options);
    end
    [rootsFound, index] = uniquetol(rootsFound, max(options.tolDelta, eps), 'DataScale', 1);
    rootStates = rootStates(index);

    function value = residualOnly(delta)
        state = evaluator(delta, fixed, options);
        if state.valid
            value = state.residual;
        else
            value = NaN;
        end
    end

end
