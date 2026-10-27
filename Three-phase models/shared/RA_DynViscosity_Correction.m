function mu_mPa_s = RA_DynViscosity_Correction(T_C, name)
    % RA_DYNVISCOSITY_CORRECTION_0 Empirical dynamic-viscosity correlations [mPa s].
    switch char(name)
        case 'Air'
            mu_mPa_s = (1.729 + 0.0049 * T_C) * 1e-2;
        case 'Water'
            mu_mPa_s = 1.3 * exp(-0.0395 * T_C) + 0.4777 * exp(-0.0061 * T_C);
        case 'Hexadecane'
            mu_mPa_s = 0.6311 * exp(239.8917 / (T_C + 106.2064));
        case 'SO-10cst'
            mu_mPa_s = 0.1057 * exp(969 / (T_C + 189.143));
        case 'SO-20cst'
            mu_mPa_s = 0.1652 * exp(1261 / (T_C + 229.1));
        otherwise
            error('Unknown fluid name: %s', string(name));
    end
end
