# Convective Epstein–Plesset Model

This folder contains the MATLAB implementation of the convective Epstein–Plesset (CEP) model used to predict the formation and temporal evolution of an interfacial gas layer around a rising oil droplet in water.

The model describes gas transfer from the surrounding water and internal oil phase into the interfacial gas layer while accounting for hydrodynamic mass transfer, hydrostatic decompression, and shear-induced gas-layer drainage.

The CEP model forms one component of the wider `gas-mediated-droplet-mobility` repository.

---

## Physical model

The model considers a spherical oil droplet rising through water while supporting an interfacial gas layer. The gas-layer volume and thickness evolve through contributions from:

- dissolved-gas transfer from the aqueous phase;
- dissolved-gas transfer from the oil phase;
- pressure variation associated with hydrostatic decompression;
- changes in the gas-layer interfacial area;
- shear-induced drainage or removal of gas from the layer.

The gas inventory is described using an ideal-gas relation, while the interfacial dissolved-gas concentrations are evaluated using Henry-law equilibrium conditions.

The water-side transport coefficient is evaluated using a Higbie-type penetration model. Oil-side gas transport is described using the corresponding internal transport formulation implemented in the model.

The model predicts an equivalent spatially averaged gas-layer thickness. It does not resolve the local angular thickness distribution around the droplet.

---

## Principal variables

The principal variables used in the CEP model include:

- `delta_a`: equivalent gas-layer thickness [m]
- `P_a`: mean gas-layer pressure [Pa]
- `V_a`: gas-layer volume [m^3]
- `n_a`: amount of gas in the interfacial layer [mol]
- `U_a`: characteristic interfacial gas-transport velocity [m/s]
- `n_dot_s`: gas transfer rate [mol/s]
- `tau_s`: characteristic shear or drainage timescale [s]
- `J_w`: water-side molar gas flux [mol/(m^2 s)]
- `J_o`: oil-side molar gas flux [mol/(m^2 s)]
- `A_w`: water-side interfacial area [m^2]
- `A_o`: oil-side interfacial area [m^2]

The precise definitions and implementation of individual variables are documented in the MATLAB functions and the accompanying research article and SI Appendix.

---

## Folder contents

- `RA_EpPlesset_Model_PNAS.m`  
  Main CEP model driver used to define the calculation and call the
  numerical solver.

- `RA_Solve_EpPlesset_PNAS.m`  
  Numerical solver for the coupled gas-layer evolution equations.

- `RA_Droplet_Correlations0.m`  
  Droplet hydrodynamic and transport correlations.

- `RA_HenryLaw.m`  
  Henry-law and gas-solubility calculations.

- `RA_HenryConcentration.m`  
  Interfacial dissolved-gas concentration calculations.

- `RA_Density_Correction_0.m`  
  Temperature-dependent density correction.

- `RA_DynViscosity_Correction_0.m`  
  Temperature-dependent dynamic-viscosity correction.

- `XSteam.m`  
  Third-party implementation of IAPWS-IF97 water and steam property calculations by Magnus Holmgren. See `XSteam_LICENSE.txt`.

---

## Running the model

1. Download or clone the complete repository.
2. Open MATLAB.
3. Set the downloaded repository as the current MATLAB folder.
4. Run the RA_Solve_EpPlesset_PNAS script.

---

## Licence

The MATLAB files beginning with `RA_` are distributed under the repository-level BSD 3-Clause License, subject to confirmation of the appropriate institutional copyright holder.

`XSteam.m` is third-party software by Magnus Holmgren and remains subject to its original BSD 3-Clause licence. See `XSteam_LICENSE.txt`.
