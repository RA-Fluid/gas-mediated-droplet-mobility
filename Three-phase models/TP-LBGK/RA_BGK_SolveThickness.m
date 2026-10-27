function result = RA_BGK_SolveThickness(fixed, options, evaluator, modelName)
    % RA_BGK_SOLVETHICKNESS Scan, bracket, solve, select, and assemble one inversion.
    [scan, brackets] = RA_BGK_ScanResidual(fixed, options, evaluator);
    if isempty(brackets)
        result = RA_BGK_AssembleResult(fixed, options, scan, [], [], 'no_bracket', modelName);
        if options.errorIfNoRoot
            error('No physical thickness root was bracketed.');
        end
        return
    end
    [rootsFound, states] = RA_BGK_SolveBrackets(brackets, fixed, options, evaluator);
    selected = RA_BGK_SelectRoot(rootsFound, states, options.rootSelection);
    result = RA_BGK_AssembleResult(fixed, options, scan, rootsFound, states, 'converged', modelName);
    result.selectedRootIndex = selected;
    result.solution = states(selected);
    result.delta_a = rootsFound(selected);
    result.success = result.solution.valid && abs(result.solution.residual) <= options.forceTolerance;
    result.diagnostics.selectedResidual_N = result.solution.residual;
    result.diagnostics.selectedRelativeResidual = result.solution.relativeResidual;
    if options.saveResult
        RA_BGK_SaveResult(result, options);
    end
end
