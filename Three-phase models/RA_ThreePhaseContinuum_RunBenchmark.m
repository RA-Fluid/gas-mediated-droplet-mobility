%% RA_ThreePhaseContinuum_RunBenchmark
clear;
clc;
close all;
caseData = struct('D_o', 343.425e-6, 'U_o', 3.159e-3, 'T_C', 22.7, ...
                'oilName', 'SO-20cst', 'rho_o', 950, 'deltaRange', [5e-9 3e-6]);
benchmark = RA_ThreePhaseContinuum_Benchmark(caseData);
if ~benchmark.allPassed
    warning('Some continuum benchmarks failed.');
end
