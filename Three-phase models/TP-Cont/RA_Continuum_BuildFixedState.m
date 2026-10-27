function fixed = RA_Continuum_BuildFixedState(D_o, U_o, T_C, oilName, options)
    % RA_CONTINUUM_BUILDFIXEDSTATE Evaluate properties independent of thickness.
    fixed = struct();
    fixed.inputs = struct('D_outer_m', D_o, 'U_d_m_s', U_o, ...
                          'T_C', T_C, 'oilName', string(oilName));
    fixed.T_K = T_C + 273.15;
    fixed.R_1 = D_o / 2;
    fixed.A_projected = pi * fixed.R_1^2;
    fixed.Psat_w = XSteam('psat_T', T_C) * 1e5;
    if options.useWaterVapourCorrection
        fixed.P_waterProperty = options.P_atm - fixed.Psat_w;
    else
        fixed.P_waterProperty = options.P_atm;
    end
    if fixed.P_waterProperty <= 0
        error('Invalid water-property pressure.');
    end
    fixed.rho_w = XSteam('rho_pT', fixed.P_waterProperty * 1e-5, T_C);
    fixed.mu_w = XSteam('my_pT', fixed.P_waterProperty * 1e-5, T_C);
    fixed.gamma_aw = XSteam('st_t', T_C);
    fixed.mu_o = RA_GetOilViscosity(T_C, oilName, options);
    fixed.rho_o = RA_GetOilDensity(T_C, oilName, options);
    fixed.mu_a = RA_GetAirViscosity(T_C, options);
    fixed.P_a = options.P_atm + fixed.rho_w * options.g * options.H + ...
        2 * fixed.gamma_aw / fixed.R_1;
    fixed.rho_a = fixed.P_a / (options.R_air * fixed.T_K);
    fixed.lambda_mfp = (fixed.mu_a / fixed.rho_a) * ...
        sqrt(pi / (2 * options.R_air * fixed.T_K));
    fixed.Re = fixed.rho_w * U_o * (2 * fixed.R_1) / fixed.mu_w;
    fixed.mu_ratio = fixed.mu_o / fixed.mu_w;
end
