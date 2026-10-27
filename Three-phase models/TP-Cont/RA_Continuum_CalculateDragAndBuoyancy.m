function state = RA_Continuum_CalculateDragAndBuoyancy(state, fixed, options)
    % RA_CONTINUUM_CALCULATEDRAGANDBUOYANCY Evaluate equivalent drag and force residual.
    state.D_U = 6 * pi * fixed.mu_w * fixed.R_1;
    state.D_s = -4 * pi * fixed.mu_w * fixed.R_1;
    state.FD_Stokes = state.D_U * fixed.inputs.U_d_m_s + state.D_s * state.A_w;
    state.f_eq = state.FD_Stokes / (2 * pi * fixed.mu_w * fixed.R_1 * fixed.inputs.U_d_m_s);
    den = 3 * fixed.mu_ratio * (2 - state.f_eq);
    if abs(den) <= 100 * eps(max(1, abs(state.f_eq)))
        state.beta_eq = Inf;
    else
        state.beta_eq = ((1 + fixed.mu_ratio) * state.f_eq - ...
                         (2 + 3 * fixed.mu_ratio)) / den;
    end
    state.b_eq = state.beta_eq * fixed.R_1;
    state.cD_Feng = (8 * state.f_eq / fixed.Re) * (1 + (state.f_eq / 16) * fixed.Re + ...
                                             (state.f_eq^2 / 160) * fixed.Re^2 * log(fixed.Re / 2));
    state.FD_Feng = 0.5 * fixed.rho_w * fixed.inputs.U_d_m_s^2 * ...
        fixed.A_projected * state.cD_Feng;
    state.FB_full = fixed.rho_w * options.g * state.V_total - ...
        fixed.rho_o * options.g * state.V_o - fixed.rho_a * options.g * state.V_a;
    state.FB_oilOnly = (fixed.rho_w - fixed.rho_o) * options.g * state.V_o;
    if strcmpi(options.buoyancyModel, 'full')
        state.FB_selected = state.FB_full;
    else
        state.FB_selected = state.FB_oilOnly;
    end
    state.cD_observed = 2 * state.FB_selected / (fixed.rho_w * ...
                                             fixed.inputs.U_d_m_s^2 * fixed.A_projected);
    state.residual = state.FD_Feng - state.FB_selected;
    state.relativeResidual = state.residual / max(abs(state.FB_selected), realmin);
end
