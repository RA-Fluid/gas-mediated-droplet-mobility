function Mol = RA_HenryConcentration(P_abs, Temp)

%--------------------------------------------------------------------------
% Main inputs
%--------------------------------------------------------------------------

% volume fractions of nitrogen (N2), oxygen (O2), Argon (Ar), and carbon dioxide (CO2) in air
volf = [0.78103, 0.20940, 0.00917, 0.00040] ; % []

% molar masses
%molm = [28.01, 32.00, 39.948, 44.01] ; % [gr mol^-1]
% molar mass of water (H2O)
molm_f = 18.01528 ; % [gr mol^-1]

% gas constant
%R = 8.3145 ; % [m^3 Pa K^-1 mol^-1]

%--------------------------------------------------------------------------
% Calculate the water saturation pressure at the given temp
%--------------------------------------------------------------------------

NT = numel(Temp) ;
NP = numel(P_abs) ;

P_abs = P_abs * 1e-3 ; % [kPa]

% water saturation pressure at the selected temperatures
Pw_sat = zeros(1, NT) ;
Pt = zeros(NP, NT) ;
for i = 1 : NT
    Pw_sat(i) = XSteam('psat_T', Temp(i)) * 100 ; % [kPa] water

    for j = 1 : NP
        Pt(j,i) = P_abs(j) - Pw_sat(i) ; % total pressure corrected by the water sta pressure at the desired temp [kPa]
    end
end

%--------------------------------------------------------------------------
% Calculate the density of the solvent (water here)
%--------------------------------------------------------------------------
densf = zeros(NP, NT) ;
for i = 1 : NP
    for j = 1 : NT
    densf(i,j) = XSteam('rho_pT', Pt(i)*0.01, Temp(j)) ; % [kg m^-3] water
    end
end

%--------------------------------------------------------------------------
% Calculate the fractions
%--------------------------------------------------------------------------
% partial pressures of nitrogen (N2), oxygen (O2), Argon (Ar) and CO2 in air

Pg_N2 = volf(1) .* Pt ; % [kPa]
Pg_O2 = volf(2) .* Pt ; % [kPa]
Pg_Ar = volf(3) .* Pt ; % [kPa]
Pg_CO2 = volf(4) .* Pt ; % [kPa]

%--------------------------------------------------------------------------
% Calculate Henry's volatility constant, Hv_px
%--------------------------------------------------------------------------
Hv_px_Pa = zeros(4, NT) ;
for i = 1: numel(Temp)
    Hv_px_Pa(1,i) = RA_HenryLaw('nitrogen', Temp(i)+273.15) ; % [Pa]
    Hv_px_Pa(2,i) = RA_HenryLaw('oxygen', Temp(i)+273.15) ; % [Pa]
    Hv_px_Pa(3,i) = RA_HenryLaw('Argon', Temp(i)+273.15) ; % [Pa]
    Hv_px_Pa(4,i) = RA_HenryLaw('carbon dioxide', Temp(i)+273.15) ; % [Pa]
end

%--------------------------------------------------------------------------
% Calculate Henry's solubility constant, Hs_cp
%--------------------------------------------------------------------------

Hs_cp_N2 = zeros(size(densf)) ;
Hs_cp_O2 = zeros(size(densf)) ;
Hs_cp_Ar = zeros(size(densf)) ;
Hs_cp_CO2 = zeros(size(densf)) ;

for i = 1 : NT
    Hs_cp_N2(:,i) = (densf(:,i) ./ (molm_f*1e-3)) / Hv_px_Pa(1,i) ; % [mol m^-3 Pa^-1]
    Hs_cp_O2(:,i) = (densf(:,i) ./ (molm_f*1e-3)) / Hv_px_Pa(2,i) ; % [mol m^-3 Pa^-1]
    Hs_cp_Ar(:,i) = (densf(:,i) ./ (molm_f*1e-3)) / Hv_px_Pa(3,i) ; % [mol m^-3 Pa^-1]
    Hs_cp_CO2(:,i) = (densf(:,i) ./ (molm_f*1e-3)) / Hv_px_Pa(4,i) ; % [mol m^-3 Pa^-1]
end

%--------------------------------------------------------------------------
% Calculate the molar fractions
%--------------------------------------------------------------------------

xg_N2 = zeros(size(densf)) ;
xg_O2 = zeros(size(densf)) ;
xg_Ar = zeros(size(densf)) ;
xg_CO2 = zeros(size(densf)) ;
for i = 1 : NP
    for j = 1 : NT
        xg_N2(i,j) = (Pg_N2(i,j)*1e3) / Hv_px_Pa(1,j) ; % []
        xg_O2(i,j) = (Pg_O2(i,j)*1e3) / Hv_px_Pa(2,j) ; % []
        xg_Ar(i,j) = (Pg_Ar(i,j)*1e3) / Hv_px_Pa(3,j) ; % []
        xg_CO2(i,j) = (Pg_CO2(i,j)*1e3) / Hv_px_Pa(4,j) ; % []
    end
end

%--------------------------------------------------------------------------
% Calculate the molarity
%--------------------------------------------------------------------------
Molg_N2 = zeros(size(densf)) ;
Molg_O2 = zeros(size(densf)) ;
Molg_Ar = zeros(size(densf)) ;
Molg_CO2 = zeros(size(densf)) ;
for j = 1 : NT
    for i = 1 : NP
        Molg_N2(i,j) = (Pg_N2(i,j)*1e3) .* Hs_cp_N2(i,j) ; % [mol m^-3]
        Molg_O2(i,j) = (Pg_O2(i,j)*1e3) .* Hs_cp_O2(i,j) ; % [mol m^-3]
        Molg_Ar(i,j) = (Pg_Ar(i,j)*1e3) .* Hs_cp_Ar(i,j) ; % [mol m^-3]
        Molg_CO2(i,j) = (Pg_CO2(i,j)*1e3) .* Hs_cp_CO2(i,j) ; % [mol m^-3]
    end
end

% molar concentration of air
Molg_Air = Molg_N2 + Molg_O2 + Molg_Ar + Molg_CO2 ;

Mol.Air = Molg_Air ;
Mol.N2 = Molg_N2 ;
Mol.O2 = Molg_O2 ;
Mol.Ar = Molg_Ar ;
Mol.CO2 = Molg_CO2 ;
