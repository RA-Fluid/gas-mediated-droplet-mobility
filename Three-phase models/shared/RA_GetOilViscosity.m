function value = RA_GetOilViscosity(T_C, oilName, options)
    if ~isempty(options.oilViscosityFcn)
        value = options.oilViscosityFcn(T_C, oilName);
    else
        value = RA_DynViscosity_Correction(T_C, oilName) * 1e-3;
    end
end
