%% RA_ThreePhaseBGK_CompareClosures.m
% Compare the closed-shell and pressure-gradient-free L-BGK closures for
% one example case

%% initilaise

clear;
clc;
close all;

%% Representative measured case

D_o = 256.6608e-6; % 200.0936e-6; %343.425e-6; % measured diameter [m]
U_o = 2.1237e-3; % 1.4719e-3; %3.159e-3; % measured rise velocity [m/s]
T_C = 22.5; % 22.7; % temperature [C]
oilName = 'SO-20cst'; % name of the oil

%% Settings common to both closures

base = struct();
base.deltaMin = 5e-9; % minimum allowable air layer thicknkess [m]
base.deltaMax = 3e-6; % maximum allowable air layer thicknkess [m]
base.nScan = 350; % number of scans to find the root
base.buoyancyModel = 'full'; % consider the full bouyancy mode always
base.saveResult = false; % do not save the results
base.errorIfNoRoot = false; % create an error and stop if there is no root

alpha_oa = 1.0; % tangential momentum accomodation coefficient; 1 is fully diffusive
alpha_aw = 1.0; % tangential momentum accomodation coefficient; 1 is fully diffusive
nVelocity = 40; % number of physicists' Hermite polynomial roots
nSpace = 101; % number of divisions in the normal direction for velocity calculation
nKnTable = 90; % number of Kn values to calculate and make the look-up table
KnMin = 0.008; % minimum Knudson number
KnMax = 2.0; % maximum Knudson number

base.kinetic = struct('alpha_oa', alpha_oa, 'alpha_aw', alpha_aw, ...
                      'nVelocity', nVelocity, 'nSpace', nSpace, 'nKnTable', nKnTable, ...
                      'KnMin', KnMin, 'KnMax', KnMax, 'interpolation', 'pchip', ...
                      'extrapolate', false, 'verbose', false);

%%  Solve both closures

% closed shell --> non-zero pressure gradient + zero net flux
closedOptions = base;
closedOptions.gasClosure = 'closed-shell';
closed = RA_ThreePhaseBGKModel(D_o, U_o, T_C, oilName, closedOptions);

% open shell --> zero pressure-gradient (Couette flow) + non-zero net flux
openOptions = base;
openOptions.gasClosure = 'zero-pressure-gradient';
openFilm = RA_ThreePhaseBGKModel(D_o, U_o, T_C, oilName, openOptions);

assert(closed.success, 'Closed-shell inversion failed.');
assert(openFilm.success, 'G = 0 inversion failed.');

%% 4. Assemble an easily readable comparison table
closure = ["Closed shell"; "G = 0 open film"];
states = {closed.solution; openFilm.solution};
comparison = table(closure, ...
                   cellfun(@(s)1e9 * s.delta_a, states), ...
                   cellfun(@(s)s.Kn, states), ...
                   cellfun(@(s)1e6 * s.b_eq, states), ...
                   cellfun(@(s)1e3 * s.A_o, states), ...
                   cellfun(@(s)1e3 * s.A_w, states), ...
                   cellfun(@(s)s.G_amplitude, states), ...
                   cellfun(@(s)s.q_relative_check, states), ...
                   cellfun(@(s)1e3 * s.meanRelativeGasVelocity, states), ...
                   cellfun(@(s)s.tau_oa_amplitude, states), ...
                   cellfun(@(s)s.tau_aw_amplitude, states), ...
                   cellfun(@(s)s.throughflowTimescale, states), ...
                   cellfun(@(s)s.relativeResidual, states), ...
                   'VariableNames', {'Closure', 'delta_nm', 'Kn', 'b_eq_um', 'Ao_mm_s', ...
                                     'Aw_mm_s', 'G_Pa_m', 'qRel_m2_s', 'uRel_mm_s', 'tau_oa_Pa', ...
                                     'tau_aw_Pa', 'throughflowTime_s', 'relativeForceResidual'});

disp(comparison);
save('RA_BGK_ClosureComparison.mat', 'comparison', 'closed', 'openFilm', '-v7.3');
