function [dydt, J_ci, J_di, q_shed] = RA_EpPlesset_Model_PNAS(t, y, p)

% y(1): Gas Film Thickness (delta)
% y(2): Internal Oil Concentration (C_oil)
% y(3): Bulk Water Concentration (C_w_infty)

delta     = y(1);
C_oil     = max(y(2), 0);
C_w_infty = max(y(3), 0);
t_eff     = max(t, 1e-6);

% HYDRODYNAMICS & ABSOLUTE PRESSURE
h_t = max(p.TotalH - (p.U * t), 0);
% h_t = p.VisH ;

P_air_total = p.Patm + (p.rho_w * p.g * h_t) + (2 * p.gamma_aw / (p.Rd + delta));
Mol = RA_HenryConcentration(P_air_total, p.T) ;

% THERMODYNAMICS
C_ci = Mol.Air ;
C_di = Mol.Air * p.K ;

% COMPOUND DROPLET KINEMATICS
% The gas film separates the water and oil, completely changing the boundary conditions.

% Outer boundary: Water shearing against Air (Free-slip limit)
% lambda_air = 1.8e-5 / p.mu_w is practically zero.
lambda_ac = p.mu_air / p.mu_w;
U_outer = p.U / (2 * (1 + lambda_ac)) ;

% Inner boundary: Air shearing against Oil
% Driven by internal circulation of the viscous core
lambda_da = p.mu_d / p.mu_air;
U_inner = p.U / (2 * (1 + lambda_da));

% The average advection velocity of the gas film
a = 1 ; % no U_inner with a = 1
U_gas_avg = a * U_outer + (1-a) * U_inner ;

tau_drain = p.Dd / U_gas_avg ; % [s]

% Volume conservation: The rate at which the gas film is swept to the wake
q_shed = delta / tau_drain ; % [m/s]

% EFFECTIVE MASS TRANSFER
% Because the water is slipping against air, the external boundary layer
% operates near the mobile Higbie limit, regardless of the oil's viscosity.
Pe_w = (p.U * p.Dd) / p.D_w;
kw_dynamic = (1.128 * sqrt(Pe_w) * p.D_w) / p.Dd;

ko_dynamic = (2 * p.D_o / p.Dd) + sqrt(p.D_o / (pi * t_eff)) ; %(p.U / (1 + lambda_da)) * 0.00375 + 

% FLUX CALCULATIONS
J_ci = kw_dynamic * (C_w_infty - C_ci);
J_di = ko_dynamic * (C_oil - C_di);
J_net = J_ci + J_di;

% ODE SYSTEM
d_delta_dt = ((p.rho_w * p.g * p.U) / P_air_total) * delta + ...
    ((p.R * (p.T+273.15)) / P_air_total) * J_net - q_shed ;

d_Coil_dt = -J_di * (3 / p.Rd);
d_Cw_infty_dt = -(4 * pi * p.Rd^2 * J_ci) / p.V_eff;

dydt = [d_delta_dt; d_Coil_dt; d_Cw_infty_dt];

end