function rho = RA_Density_Correction(T_C, name)
    % RA_DENSITY_CORRECTION_0 Empirical density correlations [kg/m^3].
    dT = T_C - 20;
    switch char(name)
        case 'Air'
            rho = 1.293 - 0.00426 * T_C;
        case 'Water'
            rho = 998.2 - 49e-4 * dT^2 - 0.2021 * dT - 0.0168;
        case 'SO-10cst'
            rho = 938.77 - 0.9216 * dT;
        case 'SO-20cst'
            rho = 952.93 - 0.9103 * dT;
        case 'Hexadecane'
            rho = 773.85 - 0.6867 * dT;
        otherwise
            error('Unknown fluid name: %s', string(name));
    end
end
