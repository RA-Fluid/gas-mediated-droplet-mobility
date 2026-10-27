function result = RA_ThreePhaseContinuumModel(D_o, U_o, T_C, oilName, userOptions)
    % RA_THREEPHASECONTINUUMMODEL Invert kinematics using the continuum gas film.
    arguments
        D_o (1, 1) double {mustBePositive}
        U_o (1, 1) double {mustBePositive}
        T_C (1, 1) double
        oilName
        userOptions struct = struct()
    end
    options = RA_BGK_MergeStruct(RA_BGK_DefaultOptions(D_o), userOptions);
    RA_BGK_ValidateOptions(options, D_o);
    fixed = RA_BGK_BuildFixedState(D_o, U_o, T_C, oilName, options);
    result = RA_BGK_SolveThickness(fixed, options, @RA_Continuum_EvaluateThickness, ...
                                   'Three-phase Maxwell-slip continuum model');
end
