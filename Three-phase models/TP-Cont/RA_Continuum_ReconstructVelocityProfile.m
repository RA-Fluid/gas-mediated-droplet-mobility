function profile = RA_Continuum_ReconstructVelocityProfile(result, nPoints)
    % RA_CONTINUUM_RECONSTRUCTVELOCITYPROFILE Exact continuum velocity profile.
    arguments
        result struct
        nPoints (1, 1) double {mustBeInteger, mustBePositive} = 201
    end
    if ~result.success
        error('The supplied solution is unsuccessful.');
    end
    s = result.solution;
    f = result.fixed;
    g = s.gas;
    eta = linspace(0, 1, nPoints).';
    y = s.delta_a * eta;
    C1C = s.DeltaA / g.L;
    C1P = -(s.G_amplitude / f.mu_a) * g.A_p / g.L;
    uC = g.l_oa * C1C + C1C * y;
    uP = g.l_oa * C1P + C1P * y + s.G_amplitude * y.^2 / (2 * f.mu_a);
    uR = uC + uP;
    uG = s.A_o + uR;
    profile = struct('eta', eta, 'y_m', y, 'y_nm', 1e9 * y, ...
                   'uCouetteRelative_m_s', uC, 'uPressure_m_s', uP, ...
                   'uRelative_m_s', uR, 'uGas_m_s', uG, ...
                   'qCouette_m2_s', trapz(y, uC), 'qPressure_m2_s', trapz(y, uP), ...
                   'qTotal_m2_s', trapz(y, uR), 'qModel_m2_s', s.q_relative, ...
                   'A_o_m_s', s.A_o, 'A_w_m_s', s.A_w, 'G_Pa_m', s.G_amplitude, ...
                   'delta_a_m', s.delta_a, 'gasClosure', result.model.gasClosureKey);
end
