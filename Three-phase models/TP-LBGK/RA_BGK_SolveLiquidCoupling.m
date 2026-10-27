function [state, ok] = RA_BGK_SolveLiquidCoupling(state, k, fixed, options)
    % RA_BGK_SOLVELIQUIDCOUPLING Solve the 2-by-2 oil/water traction system.
    ok = false;
    state.Z_o = 3 * fixed.mu_o / state.R_0;
    state.Z_w_U = 3 * fixed.mu_w / (2 * fixed.R_1);
    state.Z_w_s = -3 * fixed.mu_w / fixed.R_1;
    M = [state.Z_o - state.T_oa_Ao, -state.T_oa_Aw; ...
         -state.T_aw_Ao, state.Z_w_s - state.T_aw_Aw];
    rhs = [0; -state.Z_w_U * fixed.inputs.U_d_m_s];
    state.couplingMatrix = M;
    state.couplingRHS = rhs;
    state.matrixConditionNumber = cond(M);
    if ~all(isfinite(M), 'all') || rcond(M) <= 0 || ...
            state.matrixConditionNumber > options.maxConditionNumber
        state.message = 'Ill-conditioned kinetic liquid-coupling matrix.';
        return
    end
    A = M \ rhs;
    state.A_o = A(1);
    state.A_w = A(2);
    state.DeltaA = A(2) - A(1);
    state.G_amplitude = state.G_Ao * state.A_o + state.G_Aw * state.A_w;
    state.q_relative_check = k.qU0 * state.A_o + k.qU1 * state.A_w + k.qG * state.G_amplitude;
    state.meanRelativeGasVelocity = state.q_relative_check / state.delta_a;
    if abs(state.q_relative_check) > realmin
        state.throughflowTimescale = fixed.R_1 * state.delta_a / abs(state.q_relative_check);
    else
        state.throughflowTimescale = Inf;
    end
    state.tau_oa_amplitude = k.tau0U0 * state.A_o + k.tau0U1 * state.A_w + k.tau0G * state.G_amplitude;
    state.tau_aw_amplitude = k.tau1U0 * state.A_o + k.tau1U1 * state.A_w + k.tau1G * state.G_amplitude;
    state.tractionJumpCheck = state.tau_aw_amplitude - state.tau_oa_amplitude - ...
        state.delta_a * state.G_amplitude;
    ok = true;
end
