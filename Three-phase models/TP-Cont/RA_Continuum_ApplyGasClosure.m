function [state, ok] = RA_Continuum_ApplyGasClosure(state, gas, options)
    % RA_CONTINUUM_APPLYGASCLOSURE Apply zero-flux or zero-gradient closure.
    ok = false;
    state.gasClosure = string(options.gasClosure);
    switch lower(options.gasClosure)
        case 'closed-shell'
            if ~isfinite(gas.M_P) || gas.M_P <= 0
                state.message = 'Nonpositive pressure mobility.';
                return
            end
            state.G_over_DeltaA = gas.M_C / gas.M_P;
            Z_oa = gas.Z_C + gas.H_oa_P * state.G_over_DeltaA;
            Z_aw = gas.Z_C + gas.H_aw_P * state.G_over_DeltaA;
        case 'zero-pressure-gradient'
            state.G_over_DeltaA = 0;
            Z_oa = gas.Z_C;
            Z_aw = gas.Z_C;
        otherwise
            state.message = 'Unknown gas closure.';
            return
    end
    state.Z_oa_effective = Z_oa;
    state.Z_aw_effective = Z_aw;
    state.T_oa_Ao = -Z_oa;
    state.T_oa_Aw = Z_oa;
    state.T_aw_Ao = -Z_aw;
    state.T_aw_Aw = Z_aw;
    ok = true;
end
