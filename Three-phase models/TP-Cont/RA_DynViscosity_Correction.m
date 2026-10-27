function mu = RA_DynViscosity_Correction_0(T, name)
    % RA_DYNVISCOSITY_CORRECTION_0 Dynamic viscosity [mPa s].
    switch char(name)
        case 'Air'
            mu = (1.729 + 0.0049 * T) * 1e-2;
        case 'Water'
            mu = 1.3 * exp(-0.0395 * T) + 0.4777 * exp(-0.0061 * T);
        case 'Hexadecane'
            mu = 0.6311 * exp(239.8917 / (T + 106.2064));
        case 'SO-10cst'
            mu = 0.1057 * exp(969 / (T + 189.143));
        case 'SO-20cst'
            mu = 0.1652 * exp(1261 / (T + 229.1));
        otherwise
            error('Unknown fluid: %s', string(name));
    end
end
