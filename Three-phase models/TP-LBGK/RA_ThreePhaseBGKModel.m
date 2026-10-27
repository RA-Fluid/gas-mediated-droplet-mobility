function result = RA_ThreePhaseBGKModel(D_o, U_o, T_C, oilName, userOptions)
    % Invert droplet kinematics using a three-phase L-BGK gas film.
    %
    % D_o       measured outer air-water diameter [m]
    % U_o       measured terminal rise velocity [m/s]
    % T_C       temperature [degC]
    % oilName   'SO-20cst', 'SO-10cst', or 'Hexadecane'
    % userOptions optional structure merged recursively with defaults

    arguments
        D_o (1, 1) double {mustBePositive}
        U_o (1, 1) double {mustBePositive}
        T_C (1, 1) double
        oilName
        userOptions struct = struct()
    end

    % recall default options and merge them with user-defined options if there
    % is any the user missed to define
    options = RA_BGK_MergeStruct(RA_BGK_DefaultOptions(D_o), userOptions);

    % validate the options
    RA_BGK_ValidateOptions(options, D_o);

    % calculate properties that are independent of air layer thickness
    fixed = RA_BGK_BuildFixedState(D_o, U_o, T_C, oilName, options);

    % Build a table of outputs for the defined range of Kn if not already
    % calculated
    fixed.kineticTable = RA_BGK_LoadOrBuildResponseTable(fixed, options.kinetic);

    result = RA_BGK_SolveThickness(fixed, options, @RA_BGK_EvaluateThickness, ...
                                   'Three-phase model with linearised BGK gas film');
end
