function [rootsFound, rootStates] = RA_Continuum_SolveBrackets(brackets, fixed, options)
    % RA_CONTINUUM_SOLVEBRACKETS Solve all bracketed roots.
    count = size(brackets, 1);
    rootsFound = NaN(count, 1);
    rootStates = repmat(RA_Continuum_EmptyState(), count, 1);
    for k = 1:count
        rootsFound(k) = fzero(@residualOnly, brackets(k, :), ...
                            optimset('TolX', options.tolDelta, 'Display', 'off'));
        rootStates(k) = RA_Continuum_EvaluateThickness(rootsFound(k), fixed, options);
    end
    [rootsFound, idx] = uniquetol(rootsFound, max(options.tolDelta, eps), 'DataScale', 1);
    rootStates = rootStates(idx);

    function v = residualOnly(d)
        s = RA_Continuum_EvaluateThickness(d, fixed, options);
        if s.valid
            v = s.residual;
        else
            v = NaN;
        end
    end

end
