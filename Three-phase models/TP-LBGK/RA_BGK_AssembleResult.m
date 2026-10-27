function result = RA_BGK_AssembleResult(fixed, options, scan, roots, states, status, modelName)
    % RA_BGK_ASSEMBLERESULT Construct the public result structure.
    if strcmpi(options.gasClosure, 'closed-shell')
        closure = 'Zero relative gas flux; G solved from flux constraint';
    else
        closure = 'Pressure-gradient-free film; G=0 and qRelative reported';
    end
    result = struct();
    result.model = struct('name', modelName, 'gasClosure', closure, ...
                          'gasClosureOption', options.gasClosure, 'finiteReClosure', 'Feng level-2', ...
                          'buoyancyModel', options.buoyancyModel, 'diameterConvention', 'Outer air-water diameter');
    result.status = status;
    result.success = false;
    result.inputs = fixed.inputs;
    result.options = options;
    if isfield(fixed, 'kineticTable')
        result.fixed = rmfield(fixed, 'kineticTable');
        result.kineticTable = fixed.kineticTable;
    else
        result.fixed = fixed;
        result.kineticTable = struct();
    end
    result.scan = scan;
    result.allRoots_m = roots;
    result.allRootStates = states;
    result.solution = struct();
    result.delta_a = NaN;
    result.selectedRootIndex = NaN;
    result.diagnostics = struct('numberOfValidScanPoints', nnz([scan.valid]), ...
                                'numberOfBrackets', numel(roots), ...
                                'residualMin_N', min([scan.residual], [], 'omitnan'), ...
                                'residualMax_N', max([scan.residual], [], 'omitnan'));
end
