function state = RA_BGK_AddGeometry(state, delta, fixed)
    % Add compound-droplet geometry for a trial gas thickness.
    state.delta_a = delta;
    state.R_0 = fixed.R_1 - delta;
    state.V_total = 4 * pi * fixed.R_1^3 / 3;
    state.V_o = 4 * pi * state.R_0^3 / 3;
    state.V_a = state.V_total - state.V_o;
    state.gasVolumeFraction = state.V_a / state.V_total;
    state.effectiveDensity = (fixed.rho_o * state.V_o + fixed.rho_a * state.V_a) / state.V_total;
    state.Kn = fixed.lambda_mfp / delta;
end
