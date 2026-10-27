%% RA_ThreePhaseContinuum_CompareClosures
clear;
clc;
close all;
D_o = 343.425e-6;
U_o = 3.159e-3;
T_C = 22.7;
oilName = 'SO-20cst';
base = struct('deltaMin', 5e-9, 'deltaMax', 3e-6, 'nScan', 350, ...
            'buoyancyModel', 'full', 'saveResult', false, 'errorIfNoRoot', false);
a = base;
a.gasClosure = 'closed-shell';
b = base;
b.gasClosure = 'zero-pressure-gradient';
closedResult = RA_ThreePhaseContinuumModel(D_o, U_o, T_C, oilName, a);
openResult = RA_ThreePhaseContinuumModel(D_o, U_o, T_C, oilName, b);
if ~(closedResult.success && openResult.success)
    error('A closure failed.');
end
c = closedResult.solution;
o = openResult.solution;
Closure = ["Zero relative flux"; "Zero pressure gradient"];
delta_nm = 1e9 * [closedResult.delta_a; openResult.delta_a];
Kn = [c.Kn; o.Kn];
b_eq_um = 1e6 * [c.b_eq; o.b_eq];
Ao_mm_s = 1e3 * [c.A_o; o.A_o];
Aw_mm_s = 1e3 * [c.A_w; o.A_w];
G_Pa_m = [c.G_amplitude; o.G_amplitude];
qRel_m2_s = [c.q_relative; o.q_relative];
uRel_mm_s = 1e3 * [c.u_relative_mean; o.u_relative_mean];
tau_oa_Pa = [c.tau_oa_amplitude; o.tau_oa_amplitude];
tau_aw_Pa = [c.tau_aw_amplitude; o.tau_aw_amplitude];
advectionTime_R_s = [c.advectionTime_R; o.advectionTime_R];
relativeForceResidual = [c.relativeResidual; o.relativeResidual];
comparison = table(Closure, delta_nm, Kn, b_eq_um, Ao_mm_s, Aw_mm_s, G_Pa_m, ...
                 qRel_m2_s, uRel_mm_s, tau_oa_Pa, tau_aw_Pa, advectionTime_R_s, relativeForceResidual);
disp(comparison);
save('RA_Continuum_ClosureComparison.mat', 'comparison', 'closedResult', 'openResult');
RA_Continuum_PlotVelocityProfiles(closedResult, openResult, ...
                                  struct('saveFigure', true));
