function gas = RA_Continuum_BuildGasResponse(delta, fixed, options)
    % RA_CONTINUUM_BUILDGASRESPONSE Build Maxwell-slip Couette-Poiseuille laws.
    %
    % q_rel  = M_C*(A_w-A_o) - M_P*G
    % tau_oa = Z_C*(A_w-A_o) + H_oa_P*G
    % tau_aw = Z_C*(A_w-A_o) + H_aw_P*G
    l_oa = ((2 - options.sigma_oa) / options.sigma_oa) * fixed.lambda_mfp;
    l_aw = ((2 - options.sigma_aw) / options.sigma_aw) * fixed.lambda_mfp;
    L = delta + l_oa + l_aw;
    A_c = l_oa * delta + 0.5 * delta^2;
    A_p = 0.5 * delta^2 + l_aw * delta;
    M_C = A_c / L;
    M_P = (A_c * A_p / L - delta^3 / 6) / fixed.mu_a;
    Z_C = fixed.mu_a / L;
    H_oa_P = -A_p / L;
    H_aw_P = delta - A_p / L;
    gas = struct('l_oa', l_oa, 'l_aw', l_aw, 'L', L, 'A_c', A_c, 'A_p', A_p, ...
                 'M_C', M_C, 'M_P', M_P, 'Z_C', Z_C, ...
                 'H_oa_P', H_oa_P, 'H_aw_P', H_aw_P);
end
