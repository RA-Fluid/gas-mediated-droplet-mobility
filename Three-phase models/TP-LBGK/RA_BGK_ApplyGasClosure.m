function [state, ok] = RA_BGK_ApplyGasClosure(state, k, options)
    % RA_BGK_APPLYGASCLOSURE Reduce the kinetic response for the selected gas closure.
    ok = false;
    switch lower(options.gasClosure)
        case 'closed-shell'
            if ~isfinite(k.qG) || abs(k.qG) < realmin
                state.message = 'Degenerate kinetic pressure response.';
                return
            end
            state.G_Ao = -k.qU0 / k.qG;
            state.G_Aw = -k.qU1 / k.qG;
        case 'zero-pressure-gradient'
            state.G_Ao = 0;
            state.G_Aw = 0;
        otherwise
            state.message = sprintf('Unknown gasClosure: %s', options.gasClosure);
            return
    end
    state.T_oa_Ao = k.tau0U0 + k.tau0G * state.G_Ao;
    state.T_oa_Aw = k.tau0U1 + k.tau0G * state.G_Aw;
    state.T_aw_Ao = k.tau1U0 + k.tau1G * state.G_Ao;
    state.T_aw_Aw = k.tau1U1 + k.tau1G * state.G_Aw;
    ok = true;
end
