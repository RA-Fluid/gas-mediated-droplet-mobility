function state = RA_Continuum_EvaluateThickness(delta, fixed, options)
    % RA_CONTINUUM_EVALUATETHICKNESS Continuum Maxwell-slip calculation at one thickness.
    state = RA_BGK_EmptyState();
    state.gasClosure = string(options.gasClosure);
    if ~(isfinite(delta) && delta > 0 && delta < fixed.R_1)
        state.message = 'Nonphysical thickness.';
        return
    end
    state = RA_BGK_AddGeometry(state, delta, fixed);
    l0 = ((2 - options.sigma_oa) / options.sigma_oa) * fixed.lambda_mfp;
    l1 = ((2 - options.sigma_aw) / options.sigma_aw) * fixed.lambda_mfp;
    L = delta + l0 + l1;
    Ac = l0 * delta + 0.5 * delta^2;
    Ap = 0.5 * delta^2 + l1 * delta;
    MC = Ac / L;
    MP = (Ac * Ap / L - delta^3 / 6) / fixed.mu_a;
    if strcmpi(options.gasClosure, 'closed-shell')
        if ~isfinite(MP) || MP <= 0
            state.message = 'Nonpositive pressure mobility.';
            return
        end
        Gover = MC / MP;
    else
        Gover = 0;
    end
    Zc = fixed.mu_a / L;
    Hoa = -Ap / L;
    Haw = delta - Ap / L;
    state.G_Ao = -Gover;
    state.G_Aw = Gover;
    state.T_oa_Ao = -(Zc + Hoa * Gover);
    state.T_oa_Aw = Zc + Hoa * Gover;
    state.T_aw_Ao = -(Zc + Haw * Gover);
    state.T_aw_Aw = Zc + Haw * Gover;
    % Construct an operator-shaped structure for shared coupling/validation.
    k = struct('qU0', -MC, 'qU1', MC, 'qG', -MP, ...
               'tau0U0', -Zc, 'tau0U1', Zc, 'tau1U0', -Zc, 'tau1U1', Zc, ...
               'tau0G', Hoa, 'tau1G', Haw, 'converged', true);
    state.kinetic = k;
    [state, ok] = RA_BGK_SolveLiquidCoupling(state, k, fixed, options);
    if ~ok
        return
    end
    state = RA_BGK_CalculateDragAndBuoyancy(state, fixed, options);
    state = RA_BGK_SetValidityFlags(state, k, options);
end
