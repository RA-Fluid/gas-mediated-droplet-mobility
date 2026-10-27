%% RA_ThreePhaseBGK_RunBenchmark.m
% Verification tests for the flat modular three-phase L-BGK code.
clear;
clc;
close all;

%% 1. Canonical kinetic identities at Kn = 0.3
options = RA_BGK_DefaultOptions(300e-6);
kineticOptions = options.kinetic;
response = RA_BGK_SolvePlanarResponse(0.3, 1.0, 1.0, kineticOptions);

equalTranslationError = max(abs([ ...
                                 response.qU0 + response.qU1 - 1, ...
                                 response.tau0U0 + response.tau0U1, ...
                                 response.tau1U0 + response.tau1U1]));
symmetricCouetteError = abs(response.qU1 - 0.5);
directResidual = response.maxFixedPointError;

assert(equalTranslationError < 2e-4, 'Equal-wall translation identity failed.');
assert(symmetricCouetteError < 2e-4, 'Symmetric Couette mean failed.');
assert(directResidual < 1e-9, 'Direct kinetic-equation residual failed.');
fprintf('PASS: canonical kinetic identities.\n');

%% 2. Solve both gas closures for the representative case
RA_ThreePhaseBGK_CompareClosures;
fprintf('PASS: both closure inversions returned accepted roots.\n');
