function options = RA_BGK_DefaultOptions(D_o)
    % Central location for every model and numerical option.

    options = struct();
    options.g = 9.81;
    options.R_air = 287.05;
    options.P_atm = 101325;
    options.H = 0.11; % visualisation height
    options.rho_o = NaN;
    options.oilDensityFcn = [];
    options.oilViscosityFcn = [];
    options.airViscosityFcn = [];
    options.buoyancyModel = 'full';
    options.gasClosure = 'closed-shell';
    options.sigma_oa = 1.0;
    options.sigma_aw = 1.0;
    options.deltaMin = max(1e-9, 1e-6 * D_o);
    options.deltaMax = min(0.1 * D_o, 0.49 * D_o);
    options.nScan = 500;
    options.tolDelta = 1e-12;
    options.forceTolerance = 1e-12;
    options.rootSelection = 'smallest';
    options.maxConditionNumber = 1e12;
    options.saveResult = false;
    options.outputFile = 'RA_ThreePhaseBGKModel_Result.mat';
    options.errorIfNoRoot = false;
    options.useWaterVapourCorrection = true; % keep this always true by default

    options.kinetic = struct( ...
                             'alpha_oa', 1.0, 'alpha_aw', 1.0, ...
                             'nVelocity', 40, 'nSpace', 101, ...
                             'iterationTolerance', 2e-10, ...
                             'nKnTable', 90, 'KnMin', [], 'KnMax', [], ...
                             'interpolation', 'pchip', 'extrapolate', false, 'verbose', false);

end
