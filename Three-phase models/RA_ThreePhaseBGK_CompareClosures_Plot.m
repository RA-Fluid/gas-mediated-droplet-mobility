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

%%  Solve closed

options = base;
resultClosed = RA_ThreePhaseBGKModel( ...
                                     D_o, U_o, T_C, oilName, options);

plotOptions = struct();
plotOptions.visible = true;
plotOptions.saveFigure = true;
plotOptions.outputFile = ...
  'RA_BGK_ClosedShell_VelocityProfiles.png';

profilesClosed = RA_BGK_PlotVelocityProfiles( ...
                                             resultClosed, plotOptions);

%% solve open

options.gasClosure = 'zero-pressure-gradient';

resultOpen = RA_ThreePhaseBGKModel( ...
                                   D_o, U_o, T_C, oilName, options);

plotOptions.outputFile = ...
    'RA_BGK_OpenFilm_VelocityProfiles.png';

profilesOpen = RA_BGK_PlotVelocityProfiles( ...
                                           resultOpen, plotOptions);

%% compare

figure('Color', 'w');

plot( ...
     profilesClosed.eta, ...
     1e3 * profilesClosed.uRelative_m_s, ...
     'LineWidth', 1.8);

hold on;

plot( ...
     profilesOpen.eta, ...
     1e3 * profilesOpen.uRelative_m_s, ...
     'LineWidth', 1.8);

yline(0, ':', 'LineWidth', 1.0);

xlabel('\eta=y/\delta_a');
ylabel('u_a-A_o [mm s^{-1}]');

legend( ...
       'Zero relative flux', ...
       'Zero pressure gradient', ...
       'Location', 'best');

title('Effect of global gas-transport closure');
grid on;
box on;

%% Contour

plotOptions = struct();

plotOptions.displayThickness = 0.50;
plotOptions.nTheta = 181;
plotOptions.nEta = 81;

plotOptions.nThetaArrows = 17;
plotOptions.nEtaArrows = 9;

plotOptions.arrowScale = 0.20;

plotOptions.showAbsoluteVelocity = false;

plotOptions.saveFigure = true;
plotOptions.outputFile = ...
    'RA_BGK_SphericalRelativeVelocity.png';

sphericalOutput = ...
    RA_BGK_PlotSphericalVelocityReconstruction( ...
                                               profilesClosed, ...
                                               profilesOpen, ...
                                               plotOptions);
