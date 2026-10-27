function result = RA_Continuum_AssembleResult(fixed, options, scan, roots, states, status)
    % RA_CONTINUUM_ASSEMBLERESULT Construct public model output.
    result = struct();
    result.model = struct('name', 'Three-phase Maxwell-slip continuum model', ...
                        'gasClosureKey', options.gasClosure, 'finiteReClosure', 'Feng level-2', ...
                        'buoyancyModel', options.buoyancyModel, 'diameterConvention', 'Outer air-water diameter');
    result.status = status;
    result.success = false;
    result.inputs = fixed.inputs;
    result.options = options;
    result.fixed = fixed;
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
