function state = RA_BGK_EvaluateThickness(delta, fixed, options)
    % RA_BGK_EVALUATETHICKNESS Complete L-BGK calculation at one trial thickness.
    state = RA_BGK_EmptyState();
    state.gasClosure = string(options.gasClosure);
    if ~(isfinite(delta) && delta > 0 && delta < fixed.R_1)
        state.message = 'Nonphysical thickness.';
        return
    end
    state = RA_BGK_AddGeometry(state, delta, fixed);
    try
        k = RA_BGK_InterpolateResponse(fixed.kineticTable, state.Kn, delta, fixed, options.kinetic);
    catch ME
        state.message = ME.message;
        return
    end
    state.kinetic = k;
    [state, ok] = RA_BGK_ApplyGasClosure(state, k, options);
    if ~ok
        return
    end
    [state, ok] = RA_BGK_SolveLiquidCoupling(state, k, fixed, options);
    if ~ok
        return
    end
    state = RA_BGK_CalculateDragAndBuoyancy(state, fixed, options);
    state = RA_BGK_SetValidityFlags(state, k, options);
end
