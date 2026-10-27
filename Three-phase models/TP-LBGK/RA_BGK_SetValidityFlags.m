function state = RA_BGK_SetValidityFlags(state, k, options)
    % RA_BGK_SETVALIDITYFLAGS Apply closure-specific and common acceptance checks.
    fluxScale = max(abs(k.qU0 * state.A_o) + abs(k.qU1 * state.A_w), realmin);
    jumpScale = max([abs(state.tau_aw_amplitude), abs(state.tau_oa_amplitude), ...
                     abs(state.delta_a * state.G_amplitude), realmin]);
    state.flags = struct();
    state.flags.positiveCoreRadius = state.R_0 > 0;
    state.flags.kineticConverged = k.converged;
    state.flags.passiveFengRange = state.f_eq >= 2 && state.f_eq <= 3;
    state.flags.tractionJumpSatisfied = abs(state.tractionJumpCheck) <= 1e-2 * jumpScale;
    state.flags.finiteDrag = isfinite(state.cD_Feng) && state.cD_Feng > 0;
    if strcmpi(options.gasClosure, 'closed-shell')
        state.flags.zeroFluxSatisfied = abs(state.q_relative_check) <= 1e-8 * fluxScale;
        state.flags.zeroPressureGradientSatisfied = true;
        closureOK = state.flags.zeroFluxSatisfied;
    else
        state.flags.zeroFluxSatisfied = true;
        state.flags.zeroPressureGradientSatisfied = abs(state.G_amplitude) <= 100 * eps;
        closureOK = state.flags.zeroPressureGradientSatisfied;
    end
    state.valid = state.flags.positiveCoreRadius && state.flags.kineticConverged && ...
        state.flags.tractionJumpSatisfied && state.flags.finiteDrag && ...
        closureOK && isfinite(state.residual);
    state.message = 'ok';
end
