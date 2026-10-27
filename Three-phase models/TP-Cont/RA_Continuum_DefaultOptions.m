function options = RA_Continuum_DefaultOptions(D_o)
    % RA_CONTINUUM_DEFAULTOPTIONS Centralised defaults for the continuum model.
    options = struct();
    options.g = 9.81;
    options.R_air = 287.05;
    options.P_atm = 101325;
    options.H = 0.11;
    options.sigma_oa = 1.0;
    options.sigma_aw = 1.0;
    options.gasClosure = 'closed-shell';
    options.rho_o = NaN;
    options.oilDensityFcn = [];
    options.oilViscosityFcn = [];
    options.airViscosityFcn = [];
    options.buoyancyModel = 'full';
    options.deltaMin = max(1e-9, 1e-6 * D_o);
    options.deltaMax = min(0.1 * D_o, 0.49 * D_o);
    options.nScan = 500;
    options.tolDelta = 1e-12;
    options.forceTolerance = 1e-6;
    options.rootSelection = 'smallest';
    options.maxConditionNumber = 1e12;
    options.saveResult = false;
    options.outputFile = 'RA_ThreePhaseContinuumModel_Result.mat';
    options.errorIfNoRoot = false;
    options.useWaterVapourCorrection = true;
end
