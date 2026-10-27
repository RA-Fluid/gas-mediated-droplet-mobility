function result = RA_ThreePhaseContinuumModel(D_o, U_o, T_C, oilName, userOptions)
    % RA_THREEPHASECONTINUUMMODEL Invert gas thickness with continuum gas flow.
    %
    % Inputs: outer diameter [m], rise velocity [m/s], temperature [degC],
    % oil identifier, and an optional settings structure.
    %
    % The gas film is locally planar, uses Maxwell velocity-jump lengths,
    % and is coupled to spherical oil and water traction responses.
    arguments
        D_o (1, 1) double {mustBePositive}
        U_o (1, 1) double {mustBePositive}
        T_C (1, 1) double
        oilName
        userOptions struct = struct()
    end

    options = RA_Continuum_DefaultOptions(D_o);
    options = RA_Continuum_MergeStruct(options, userOptions);
    RA_Continuum_ValidateOptions(options, D_o);
    fixed = RA_Continuum_BuildFixedState(D_o, U_o, T_C, oilName, options);
    [scan, brackets] = RA_Continuum_ScanResidual(fixed, options);

    if isempty(brackets)
        result = RA_Continuum_AssembleResult(fixed, options, scan, [], [], 'no_bracket');
        if options.saveResult
            RA_Continuum_SaveResult(result, options);
        end
        if options.errorIfNoRoot
            error('RA_Continuum:NoBracket', 'No thickness root was bracketed.');
        end
        return
    end

    [rootsFound, rootStates] = RA_Continuum_SolveBrackets(brackets, fixed, options);
    selected = RA_Continuum_SelectRoot(rootsFound, rootStates, options.rootSelection);
    result = RA_Continuum_AssembleResult(fixed, options, scan, rootsFound, rootStates, 'converged');
    result.selectedRootIndex = selected;
    result.solution = rootStates(selected);
    result.delta_a = rootsFound(selected);
    result.success = result.solution.valid && ...
        abs(result.solution.relativeResidual) <= options.forceTolerance;
    result.diagnostics.selectedResidual_N = result.solution.residual;
    result.diagnostics.selectedRelativeResidual = result.solution.relativeResidual;
    if options.saveResult
        RA_Continuum_SaveResult(result, options);
    end
end
