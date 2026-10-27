function value = RA_GetOilDensity(T_C, oilName, options)

    if isfinite(options.rho_o)
        value = options.rho_o;
    elseif ~isempty(options.oilDensityFcn)
        value = options.oilDensityFcn (T_C, oilName);
    else
        value = RA_Density_Correction(T_C, oilName);
    end
