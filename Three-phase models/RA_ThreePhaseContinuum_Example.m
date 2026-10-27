%% RA_ThreePhaseContinuum_Example.m
% Run the three-phase continuum model for one representative observation.
%
% Select either:
%   options.gasClosure = 'closed-shell';
%   options.gasClosure = 'zero-pressure-gradient';

clear;
clc;
close all;

% Get the directory path of this main script
scriptFolder = fileparts(mfilename('fullpath'));

% Add model subfolders to the MATLAB search path
addpath(fullfile(scriptFolder, 'shared'));
addpath(fullfile(scriptFolder, 'TP-Cont'));

%% 1. Experimental input

D_o = 343.425e-6;       % Measured outer droplet diameter [m]
U_o = 3.159e-3;         % Measured terminal rise velocity [m/s]
T_C = 22.7;             % Temperature [degC]
oilName = 'SO-20cst';   % Oil identifier

%% 2. Model and numerical options

options = struct();

% Gas-transport closure:
%
% 'closed-shell'
%     Enforce zero relative tangential gas flux and solve for G.
%
% 'zero-pressure-gradient'
%     Set G = 0 and calculate the resulting relative gas flux.
options.gasClosure = 'closed-shell'; %'zero-pressure-gradient'

% Thickness-search interval.
options.deltaMin = 5e-9;       % Minimum trial thickness [m]
options.deltaMax = 3e-6;       % Maximum trial thickness [m]
options.nScan = 350;           % Number of logarithmic scan points

% Physical and output settings.
options.buoyancyModel = 'full';
options.saveResult = false;
options.errorIfNoRoot = false;

%% 3. Run the model

result = RA_ThreePhaseContinuumModel( ...
                                     D_o, ...
                                     U_o, ...
                                     T_C, ...
                                     oilName, ...
                                     options);

%% 4. Display the solver status

fprintf('\n');
fprintf('THREE-PHASE CONTINUUM MODEL\n');
fprintf('===========================\n');
fprintf('Status: %s\n', result.status);

%% 5. Print the converged result

if result.success

    state = result.solution;

    fprintf('\n');
    fprintf('Selected gas closure\n');
    fprintf('--------------------\n');
    fprintf('%s\n', result.model.gasClosureKey);

    fprintf('\n');
    fprintf('Inferred quantities\n');
    fprintf('-------------------\n');
    fprintf('Gas-layer thickness, delta_a : %.3f nm\n', ...
            1e9 * result.delta_a);

    fprintf('Gas-film Knudsen number, Kn  : %.6g\n', ...
            state.Kn);

    fprintf('Equivalent slip length, b_eq : %.3f um\n', ...
            1e6 * state.b_eq);

    fprintf('\n');
    fprintf('Interfacial velocity amplitudes\n');
    fprintf('-------------------------------\n');
    fprintf('Oil-air amplitude, A_o        : %.6g mm/s\n', ...
            1e3 * state.A_o);

    fprintf('Air-water amplitude, A_w      : %.6g mm/s\n', ...
            1e3 * state.A_w);

    fprintf('Relative amplitude, DeltaA    : %.6g mm/s\n', ...
            1e3 * state.DeltaA);

    fprintf('\n');
    fprintf('Gas transport\n');
    fprintf('-------------\n');
    fprintf('Pressure-gradient amplitude, G: %.6e Pa/m\n', ...
            state.G_amplitude);

    fprintf('Relative gas flux, q_rel      : %.6e m^2/s\n', ...
            state.q_relative);

    fprintf('Mean relative gas velocity    : %.6g mm/s\n', ...
            1e3 * state.u_relative_mean);

    fprintf('Radius-based advection time   : %.6g s\n', ...
            state.advectionTime_R);

    fprintf('\n');
    fprintf('Gas tractions\n');
    fprintf('-------------\n');
    fprintf('Oil-air traction, tau_oa      : %.6e Pa\n', ...
            state.tau_oa_amplitude);

    fprintf('Air-water traction, tau_aw    : %.6e Pa\n', ...
            state.tau_aw_amplitude);

    fprintf('Traction-jump residual        : %.6e Pa\n', ...
            state.tractionJumpCheck);

    fprintf('\n');
    fprintf('Force-balance diagnostics\n');
    fprintf('-------------------------\n');
    fprintf('Predicted drag, F_D            : %.6e N\n', ...
            state.FD_Feng);

    fprintf('Selected buoyancy, F_B         : %.6e N\n', ...
            state.FB_selected);

    fprintf('Absolute force residual        : %.6e N\n', ...
            state.residual);

    fprintf('Relative force residual        : %.6e\n', ...
            state.relativeResidual);

    fprintf('Coupling-matrix condition no.  : %.6e\n', ...
            state.matrixConditionNumber);

    fprintf('\n');
    fprintf('Closure checks\n');
    fprintf('--------------\n');
    fprintf('Closure condition satisfied    : %d\n', ...
            state.flags.closureSatisfied);

    fprintf('Traction identity satisfied    : %d\n', ...
            state.flags.tractionJumpSatisfied);

    fprintf('Finite positive drag           : %d\n', ...
            state.flags.finiteDrag);

else

    fprintf('\n');
    warning([ ...
             'No accepted gas-layer-thickness root was found. ', ...
             'Inspect result.scan and result.diagnostics.']);

    fprintf('Valid scan points : %d\n', ...
            result.diagnostics.numberOfValidScanPoints);

    fprintf('Root brackets      : %d\n', ...
            result.diagnostics.numberOfBrackets);

    fprintf('Minimum residual   : %.6e N\n', ...
            result.diagnostics.residualMin_N);

    fprintf('Maximum residual   : %.6e N\n', ...
            result.diagnostics.residualMax_N);
end
