%--------------------------------------------------------------------------
% This scripts solves the Epstein-Plesset equation developed to estimate
% an air layer thickness around a rising droplet in a liquid column, based
% on the thermodynamic saturation state of both phases
%
% Developed by Reza Azadi (2026), The University of Sydney, School of
% Chemistry
%
%--------------------------------------------------------------------------

%% initiate

clear
clc
close all

%% main inputs

% Atmospheric Pressure 
datInput.Patm = 101325; % [Pa]
% Temperature
datInput.T = 22; % [C]

% bulk water saturation level
% this extra 0.035 i sto account for the saturation level under the
% hydrostatic pressure; meaning that 1.035 is actually 1.00 in ambient
% because pressure is higher at bottom
Si_W  = 1.2 + 0.035 ;
% Oil initial saturation level
Si_O = 1.0 + 0.035 ;

datInput.Si_W = Si_W ;
datInput.Si_O = Si_O ;

% diameter range
Dd_vec = (150:10:400)' * 1e-6; % [m]

% name of the lubricants
oilNames = {'Hexadecane'; 'SO-20cst' } ;

% Selected diamater to plot the temporal evolution of air layer thickness
Dd_sel = 250 ; % [um]

%% main fixed variables

% Universal Gas Constant 
datInput.R = 8.31446; % [J/(mol.K)]
% Gravity
datInput.g = 9.81; % [m/s^2]
% air dynamic viscosity
datInput.mu_air = RA_DynViscosity_Correction_0(datInput.T, 'Air') * 1e-3 ; % [Pa.s]

% Liquid column geometry
datInput.TotalH = 0.13; % Total liquid column height (m)
datInput.VisH = 0.11; % Height of the visualization plane [m]
datInput.Vol_Tank = 0.06 * 0.06 * datInput.TotalH ; % Total water volume [mL]

%% Water properties

% saturation pressure of water at the defined temperature
datInput.Psat_w = XSteam( 'psat_T', datInput.T ) * 1e5 ; % [Pa]
% total pressure of water at the defined temperature
datInput.Ptot_w = datInput.Patm - datInput.Psat_w ; % [Pa]
% density of the water
datInput.rho_w = XSteam('rho_pT', datInput.Ptot_w*1e-5, datInput.T ) ; % [kg/m^3]
% Dynamic viscosity of water
datInput.mu_w = XSteam('my_pT', datInput.Ptot_w*1e-5, datInput.T ) ; % [Pa.s]
% Air-Water surface tension
datInput.gamma_aw = XSteam('st_t', datInput.Ptot_w*1e-5, datInput.T ) ; % [N/m]

%--------------------------------------------------------------------------
% Diffusivity of water
%--------------------------------------------------------------------------
% based on the Diaz et al's (1987) correlation [20-30 deg C]
mu_25C = XSteam('my_pT', datInput.Ptot_w*1e-5, 25 ) ; % [Pa.s]
T_Kelvin = datInput.T + 273.15 ;
D_w = 2.7e-9 * (T_Kelvin/298.15) * (mu_25C / datInput.mu_w) ; % [m^2/s]

%--------------------------------------------------------------------------
% Air solubility in water
%--------------------------------------------------------------------------
Mol = RA_HenryConcentration(datInput.Patm, datInput.T) ;
datInput.Molarity_Air = Mol.Air ; % [mol / m^3]
datInput.Cw_infty_init = Si_W * datInput.Molarity_Air; % [mol / m^3]

%% Oil properties
%--------------------------------------------------------------------------
% surface tension
%--------------------------------------------------------------------------
Hexadecane.a_0 = 27.15 ; % [mN/m]
Hexadecane.a_1 = 0.086 ; % [mN/(m.K)]

SO_PDMS.a_0 = 20.10 ; % [mN/m]
SO_PDMS.a_1 = 0.065 ; % [mN/(m.K)]

Hexadecane.sigma = Hexadecane.a_0 - Hexadecane.a_1 * (T_Kelvin - 298.15) ; % [mN/m]
SO_PDMS.sigma = SO_PDMS.a_0 - SO_PDMS.a_1 * (T_Kelvin - 298.15) ; % [mN/m]

%--------------------------------------------------------------------------
% dynamic viscosity: from measurements
%--------------------------------------------------------------------------
Hexadecane.mu = RA_DynViscosity_Correction_0(datInput.T, oilNames{1}) ; % [mPa.s]
SO_PDMS.mu = RA_DynViscosity_Correction_0(datInput.T,  oilNames{2}) ; % [mPa.s]

%--------------------------------------------------------------------------
% Diffusivity of oils
%--------------------------------------------------------------------------
Va = 29.9 ; % [cm^3/mol]
eps = (10.2/Va) - 0.791 ;
Hexadecane.D_o = (13.3e-12) * (T_Kelvin^1.47) * (Hexadecane.mu^eps) * (Va^-0.71) ; % [m^2/s]
SO_PDMS.D_o = (13.3e-12) * (T_Kelvin^1.47) * (SO_PDMS.mu^eps) * (Va^-0.71) ; % [m^2/s]

%--------------------------------------------------------------------------
% density: from measurements
%--------------------------------------------------------------------------
Hexadecane.rho = RA_Density_Correction_0(datInput.T,  oilNames{1}) ;
SO_PDMS.rho = RA_Density_Correction_0(datInput.T,  oilNames{2}) ;

%% Allocate the arrays

nD = length(Dd_vec);
nOil = size(oilNames, 1) ;

% Preallocate Data Matrices
deltaAirFinal = zeros(nOil, nD);
timeFinal = zeros(nOil, nD);
lengthFinal = zeros(nOil, nD);

%% Velocities | Schiller-Naumann

Corrs = RA_Droplet_Correlations0(Hexadecane.sigma, datInput.rho_w, datInput.mu_w, Hexadecane.rho, Hexadecane.mu*1e-3, Dd_vec) ;
Hexadecane.U_d = Corrs.SolidSphere.Schiller_Numann.U_d ; % [mm/s]

Corrs = RA_Droplet_Correlations0(SO_PDMS.sigma, datInput.rho_w, datInput.mu_w, SO_PDMS.rho, SO_PDMS.mu*1e-3, Dd_vec) ;
SO_PDMS.U_d = Corrs.SolidSphere.Schiller_Numann.U_d ; % [mm/s]

%% Run the solver in the loop

fprintf('Running diameter Sweep (100 to 500 um)...\n');

% Preparefigure
figure;
hold on
grid on
box on

subplot(1,3,1) ;
title( ['$\rm Dynamic \, Gas \, Film \, Evolution \, (D_d = 250\,', '\mu', '\rm m)$'], 'Interpreter', 'latex', 'FontSize', 16);
grid on
box on
hold on
xlabel('Rise Time (s)', 'Interpreter', 'latex', 'FontSize', 14) ; 
ylabel('Gas Film Thickness, $\delta_{\rm air}$ (nm)', 'Interpreter', 'latex', 'FontSize', 14) ;
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex'); 
ax = gca ;
ax.XScale = 'log' ;

RA_Output = struct() ;
for i = 1:nOil
    
    for j = 1:nD

        Dd = Dd_vec(j);
        Rd = Dd / 2;           
        Vol_drop = (4/3)*pi*Rd^3;
        
        if contains(oilNames{i,1}, 'Hexadecane')

            K = 4.2 ;
            U = Hexadecane.U_d(j)*1e-3 ; % [m/s]
            D_o = Hexadecane.D_o ; % [m^2/s]
            mu_d = Hexadecane.mu * 1e-3 ; % [Pa.s]
            sigma_d = Hexadecane.sigma*1e-3 ; % [N/m]
            rho_d = Hexadecane.rho ; % [kg/m3]

        elseif contains(oilNames{i,1}, 'SO')

            K = 6.0 ;
            U = SO_PDMS.U_d(j)*1e-3 ; % [m/s]
            D_o = SO_PDMS.D_o ; % [m^2/s]
            mu_d = SO_PDMS.mu * 1e-3 ; % [Pa.s]
            sigma_d = SO_PDMS.sigma*1e-3 ; % [N/m]
            rho_d = SO_PDMS.rho ; % [kg/m3]
        end
        
        % Parameters
        par = datInput;
        par.mu_d = mu_d;
        par.sigma_d = sigma_d ;
        par.rho_d = rho_d ;
        par.D_w = D_w;
        par.D_o = D_o;
        par.Dd = Dd;
        par.Rd = Rd;
        par.U = U;       
        par.K = K;
        par.Vol_drop = Vol_drop;
        par.V_eff = datInput.Vol_Tank;
        par.Co_init = Si_O * datInput.Molarity_Air * K;
        
        % Solve up to the visualisation window
        t_end = datInput.VisH / U; % [s]
        
        % set the solver options
        options = odeset('RelTol', 1e-6, 'AbsTol', [1e-12; 1e-4; 1e-4], ...
                         'NonNegative', [1, 2, 3], 'MaxStep', t_end/500);
        
        % initial conditions
        y0 = [1e-9; par.Co_init; datInput.Cw_infty_init];
        
        % Run ODE
        [t, y] = ode15s(@(t,y) RA_EpPlesset_Model_PNAS(t, y, par), [1e-4 t_end], y0, options);

        %------------------------------------------------------------------
        % calculate the instantaneous fluxes
        %------------------------------------------------------------------
        % pre-allocate arrays for the results of this specific run
        nSteps = length(t);
        J_ci_vec = zeros(nSteps, 1);
        J_di_vec = zeros(nSteps, 1);
        q_shed_vec = zeros(nSteps, 1);

        % loop through the solver output to re-calculate the internal variables
        for k = 1:nSteps
            [~, J_ci_vec(k), J_di_vec(k), q_shed_vec(k)] = RA_EpPlesset_Model_PNAS(t(k), y(k,:)', par);
        end
        
        % the air layer thickness
        delta_nm = y(:,1) * 1e9;
        
        % extract Final Thickness at the visualisation height
        delta_inf = delta_nm(end);
        deltaAirFinal(i, j) = delta_inf;
        
        % plot on Fig 1 if Dd is exactly Dd_sel um
        if abs(Dd - Dd_sel *1e-6) < 1e-9
            plot(t, delta_nm, 'LineWidth', 2.5, 'DisplayName', oilNames{i,1});
        end
        
        % calculate Time and Length to the visualisation height
        timeFinal(i, j) = t(end);
        lengthFinal(i, j) = t(end) * U * 1000; % Length in mm

        RA_Output(j,i).datInput = par ;
        RA_Output(j,i).t = t ;
        RA_Output(j,i).delta_nm = delta_nm ;
        RA_Output(j,i).J_ci = J_ci_vec ;
        RA_Output(j,i).J_di = J_di_vec ;
        RA_Output(j,i).q_shed = q_shed_vec ;
        RA_Output(j,i).delta_inf = delta_inf ;
        RA_Output(j,i).timeFinal = t(end) ;
        RA_Output(j,i).lengthFinal = t(end) * U * 1000 ;

    end
end
hold off
legend('Location','best', 'Interpreter', 'latex', 'FontSize', 12); 
drawnow

%% Final Thickness vs Diameter

subplot(1,3,2) ;
grid on
box on
hold on
colors = lines(nOil);
for i = 1:nOil
    plot(Dd_vec*1e6, deltaAirFinal(i,:), '-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
        'DisplayName', oilNames{i,1}, 'Color', colors(i,:), 'MarkerFaceColor', 'none');
end
hold off

xlabel('Droplet Diameter, $D_d$ ($\mu$m)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Air Layer Thickness, $\delta_{\rm a}$ (nm)', 'Interpreter', 'latex', 'FontSize', 14);
title('Equilibrium Gas Film Thickness vs. Droplet Size', 'Interpreter', 'latex', 'FontSize', 16);
legend('Location','northeast', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
drawnow

%% Time and length scales to visualisation height
subplot(1,3,3) ;
grid on
box on
hold on

for i = 1:nOil
    plot(Dd_vec*1e6, timeFinal(i,:), '-', 'LineWidth', 2, 'MarkerSize', 7, ...
        'DisplayName', [oilNames{i,1} ' (Time)']);
end
hold off

% Formatting for Dual Axis
ylabel('Time to reach 110 mm (s)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Droplet Diameter, $D_d$ ($\mu$m)', 'Interpreter', 'latex', 'FontSize', 14);
% title('Kinematic Scales for Film Equilibration', 'Interpreter', 'latex', 'FontSize', 16);

% Creating a custom legend to handle both axes cleanly
legend('Location','northwest', 'Interpreter', 'latex', 'FontSize', 10);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
drawnow

fprintf('Sweep Complete. Figures generated.\n');

%% save the output

% create the filename
baseName = sprintf('RA_Output_Si_W_%.2f_Si_O_%.2f.mat', Si_W-0.035, Si_O-0.035);
filname = fullfile('output', baseName);

baseName = sprintf('RA_Output_Si_W_%.2f_Si_O_%.2f.svg', Si_W-0.035, Si_O-0.035);
Imname =  fullfile('output', baseName);

% Check if the directory exists; if not, create it to avoid a save error
if ~exist('output', 'dir')
    mkdir('output');
end

% Save the structure
save(filname, 'RA_Output', '-mat');
exportgraphics(gcf, Imname, 'Resolution', 600) ;