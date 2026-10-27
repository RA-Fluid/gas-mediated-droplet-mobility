function [state, ok] = RA_Continuum_SolveLiquidCoupling(state, gas, fixed, options)
    % RA_CONTINUUM_SOLVELIQUIDCOUPLING Solve the coupled interface amplitudes.
    %
    % The oil-side and water-side continuum traction responses are coupled to
    % the reduced gas-film traction operator after applying the selected gas
    % closure. The unknowns are the oil-air and air-water tangential velocity
    % amplitudes, A_o and A_w.
    %
    % IMPORTANT
    %   MATLAB matrix left division must be written exactly as:
    %
    %       amplitudes = matrix \ rhs;
    %
    %   The previous generated file contained a corrupted line break at this
    %   statement. This replacement file corrects that generation error.

    ok = false;

    % Continuum liquid traction coefficients.
    state.Z_o   =  3 * fixed.mu_o / state.R_0;
    state.Z_w_U =  3 * fixed.mu_w / (2 * fixed.R_1);
    state.Z_w_s = -3 * fixed.mu_w / fixed.R_1;

    % Coupled traction-balance system.
    % Gas tractions after closure are written as:
    %   tau_oa = T_oa_Ao*A_o + T_oa_Aw*A_w
    %   tau_aw = T_aw_Ao*A_o + T_aw_Aw*A_w
    matrix = [ ...
              state.Z_o - state.T_oa_Ao,  -state.T_oa_Aw; ...
              -state.T_aw_Ao,              state.Z_w_s - state.T_aw_Aw];

    rhs = [ ...
           0; ...
           -state.Z_w_U * fixed.inputs.U_d_m_s];

    state.couplingMatrix = matrix;
    state.couplingRHS = rhs;
    state.matrixConditionNumber = cond(matrix);

    if ~all(isfinite(matrix), 'all') || rcond(matrix) <= 0 || ...
            state.matrixConditionNumber > options.maxConditionNumber
        state.message = 'Ill-conditioned liquid-coupling matrix.';
        return
    end

    % Solve simultaneously for A_o and A_w.
    amplitudes = matrix \ rhs;

    state.A_o = amplitudes(1);
    state.A_w = amplitudes(2);
    state.DeltaA = state.A_w - state.A_o;

    % Recover the closure-dependent pressure gradient and relative gas flux.
    state.G_amplitude = state.G_over_DeltaA * state.DeltaA;
    state.q_relative = gas.M_C * state.DeltaA - ...
        gas.M_P * state.G_amplitude;
    state.u_relative_mean = state.q_relative / state.delta_a;

    % Recover the separate gas tractions at both liquid interfaces.
    state.tau_oa_amplitude = gas.Z_C * state.DeltaA + ...
        gas.H_oa_P * state.G_amplitude;
    state.tau_aw_amplitude = gas.Z_C * state.DeltaA + ...
        gas.H_aw_P * state.G_amplitude;

    % Continuum tangential momentum identity:
    %   tau_aw - tau_oa = delta_a*G.
    state.tractionJumpCheck = state.tau_aw_amplitude - ...
        state.tau_oa_amplitude - state.delta_a * state.G_amplitude;

    % Characteristic radius-based tangential advection time. This is a local
    % transport scale, not a global gas residence or circulation time.
    if abs(state.u_relative_mean) > realmin
        state.advectionTime_R = fixed.R_1 / abs(state.u_relative_mean);
    else
        state.advectionTime_R = Inf;
    end

    ok = true;
end
