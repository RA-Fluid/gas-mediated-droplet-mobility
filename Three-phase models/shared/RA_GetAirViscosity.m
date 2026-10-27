function value = RA_GetAirViscosity(T_C, options)

    if ~isempty(options.airViscosityFcn)
        value = options.airViscosityFcn(T_C);
    else
        value = RA_DynViscosity_Correction(T_C, 'Air') * 1e-3;
    end
