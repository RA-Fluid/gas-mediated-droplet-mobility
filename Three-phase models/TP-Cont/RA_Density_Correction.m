function rho = RA_Density_Correction_0(T, name)
    % RA_DENSITY_CORRECTION_0 Density [kg/m^3].
    dT = T - 20;
    switch char(name)
        case 'Air'
            rho = 1.293 - 0.00426 * T;
        case 'Water'
            rho = 998.2 - 49e-4 * dT^2 - 0.2021 * dT - 0.0168;
        case 'SO-10cst'
            rho = 938.77 - 0.9216 * dT;
        case 'SO-20cst'
            rho = 952.93 - 0.9103 * dT;
        case 'Hexadecane'
            rho = 773.85 - 0.6867 * dT;
        otherwise
            error('Unknown fluid: %s', string(name));
    end
end
