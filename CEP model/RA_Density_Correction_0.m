function rho = RA_Density_Correction_0(T, name)

T0 = 20 ; % [c]
dT = T - T0 ; % [c]

switch name

    case 'Air'
        rho = 1.293 - 0.00426 * T ; % [kg/m3]
        
    case 'Water'
        % rho_0 = 998.2 ; % [kg/m3]
        % rho = rho_0 -49e-4*(dT^2) - 0.2021 * dT - 0.0168 ; % [kg/m3]
        rho = XSteam('rho_PT', 1.01325, T) ; % [kg/m3]

    case 'SO-10cst'
        rho_0 = 938.77 ; % [kg/m3]
        rho = rho_0 -0.9216*dT ; % [kg/m3]

    case 'SO-20cst'
        rho_0 = 952.93 ; % [kg/m3]
        rho = rho_0 -0.9103*dT ; % [kg/m3]

    case 'Hexadecane'
        rho_0 = 773.85 ; % [kg/m3]
        rho = rho_0 -0.6867*dT ; % [kg/m3]
end