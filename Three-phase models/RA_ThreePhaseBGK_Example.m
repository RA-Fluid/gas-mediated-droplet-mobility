%% RA_ThreePhaseBGK_Example.m
% Run one three-phase L-BGK inversion using either gas-film closure.
clear;
clc;
close all;

% Get the directory path of this main script
scriptFolder = fileparts(mfilename('fullpath'));

% Add model subfolders to the MATLAB search path
addpath(fullfile(scriptFolder, 'shared'));
addpath(fullfile(scriptFolder, 'TP-LBGK'));

%% 1. Measured droplet data
D_o = 343.425e-6;              % outer droplet diameter [m]
U_o = 3.159e-3;                % terminal rise velocity [m/s]
T_C = 22.7;                    % temperature [degC]
oilName = 'SO-20cst';

%% 2. Numerical and physical options
options = struct();
% options.gasClosure = 'closed-shell';
% Alternative pressure-gradient-free limit:
options.gasClosure = 'zero-pressure-gradient';
options.deltaMin = 5e-9;
options.deltaMax = 3e-6;
options.nScan = 350;
options.buoyancyModel = 'full';
options.saveResult = false;
options.errorIfNoRoot = false;
options.kinetic = struct();
options.kinetic.alpha_oa = 1.0;
options.kinetic.alpha_aw = 1.0;
options.kinetic.nVelocity = 40;
options.kinetic.nSpace = 101;
options.kinetic.nKnTable = 90;
options.kinetic.KnMin = 0.008;
options.kinetic.KnMax = 2.0;
options.kinetic.interpolation = 'pchip';
options.kinetic.extrapolate = false;
options.kinetic.verbose = true;

%% 3. Run the inversion
result = RA_ThreePhaseBGKModel(D_o, U_o, T_C, oilName, options);
disp(result.status);

%% 4. Print the selected physical root
if result.success
    s = result.solution;
    fprintf('\nClosure: %s\n', options.gasClosure);
    fprintf('delta_a = %.3f nm\n', 1e9 * s.delta_a);
    fprintf('Kn = %.6g\n', s.Kn);
    fprintf('b_eq = %.3f um\n', 1e6 * s.b_eq);
    fprintf('A_o = %.6g mm/s\n', 1e3 * s.A_o);
    fprintf('A_w = %.6g mm/s\n', 1e3 * s.A_w);
    fprintf('G = %.6e Pa/m\n', s.G_amplitude);
    fprintf('q_relative = %.6e m^2/s\n', s.q_relative_check);
    fprintf('mean relative gas velocity = %.6e mm/s\n', 1e3 * s.meanRelativeGasVelocity);
    fprintf('throughflow time = %.6e s\n', s.throughflowTimescale);
    fprintf('tau_oa = %.6e Pa\n', s.tau_oa_amplitude);
    fprintf('tau_aw = %.6e Pa\n', s.tau_aw_amplitude);
    fprintf('traction jump residual = %.3e Pa\n', s.tractionJumpCheck);
    fprintf('relative force residual = %.3e\n', s.relativeResidual);
else
    warning('No accepted root. Inspect result.scan and result.diagnostics.');
end
