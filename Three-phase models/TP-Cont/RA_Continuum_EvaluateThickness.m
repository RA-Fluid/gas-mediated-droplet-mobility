function state = RA_Continuum_EvaluateThickness(delta, fixed, options)
    % RA_CONTINUUM_EVALUATETHICKNESS Evaluate one trial gas-layer thickness.
    state = RA_Continuum_EmptyState();
    state.delta_a = delta;
    if ~(isfinite(delta) && delta > 0 && delta < fixed.R_1)
        state.message = 'Nonphysical thickness.';
        return
    end
    state.R_0 = fixed.R_1 - delta;
    state.V_total = 4 * pi * fixed.R_1^3 / 3;
    state.V_o = 4 * pi * state.R_0^3 / 3;
    state.V_a = state.V_total - state.V_o;
    state.gasVolumeFraction = state.V_a / state.V_total;
    state.effectiveDensity = (fixed.rho_o * state.V_o + fixed.rho_a * state.V_a) / state.V_total;
    state.Kn = fixed.lambda_mfp / delta;
    gas = RA_Continuum_BuildGasResponse(delta, fixed, options);
    state.gas = gas;
    [state, ok] = RA_Continuum_ApplyGasClosure(state, gas, options);
    if ~ok
        return
    end
    [state, ok] = RA_Continuum_SolveLiquidCoupling(state, gas, fixed, options);
    if ~ok
        return
    end
    state = RA_Continuum_CalculateDragAndBuoyancy(state, fixed, options);
    fluxScale = max(abs(gas.M_C * state.DeltaA), realmin);
    jumpScale = max([abs(state.tau_aw_amplitude), abs(state.tau_oa_amplitude), ...
                     abs(delta * state.G_amplitude), realmin]);
    state.flags = struct();
    state.flags.positiveCoreRadius = state.R_0 > 0;
    state.flags.positiveMobilities = gas.M_C > 0 && gas.M_P > 0;
    state.flags.finiteDrag = isfinite(state.cD_Feng) && state.cD_Feng > 0;
    state.flags.tractionJumpSatisfied = abs(state.tractionJumpCheck) <= 1e-9 * jumpScale;
    if strcmpi(options.gasClosure, 'closed-shell')
        state.flags.closureSatisfied = abs(state.q_relative) <= 1e-9 * fluxScale;
    else
        state.flags.closureSatisfied = abs(state.G_amplitude) <= 1e-14;
    end
    state.valid = state.flags.positiveCoreRadius && state.flags.positiveMobilities && ...
        state.flags.finiteDrag && isfinite(state.residual) && ...
        state.matrixConditionNumber <= options.maxConditionNumber;
    if state.valid
        state.message = 'ok';
    end
end
