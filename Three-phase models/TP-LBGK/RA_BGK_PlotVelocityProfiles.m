function profiles = RA_BGK_PlotVelocityProfiles(result, plotOptions)
    % RA_BGK_PLOTVELOCITYPROFILES Reconstruct and plot the converged gas velocity.
    %
    % PURPOSE
    %   Reconstruct the dimensional gas velocity profile from the three
    %   linearised-BGK basis solutions:
    %
    %       1. unit oil-air interface motion;
    %       2. unit air-water interface motion;
    %       3. unit pressure-driving source.
    %
    %   The function plots:
    %       - gas velocity in the adopted laboratory/droplet frame;
    %       - gas velocity relative to the oil-air interface;
    %       - Couette and pressure-driven contributions;
    %       - cumulative relative gas flux across the film.
    %
    % INPUTS
    %   result
    %       Successful output from RA_ThreePhaseBGKModel.
    %
    %   plotOptions
    %       Optional plotting settings:
    %
    %       plotOptions.visible
    %           true or false. Default: true.
    %
    %       plotOptions.saveFigure
    %           true or false. Default: false.
    %
    %       plotOptions.outputFile
    %           Figure filename. Default:
    %           'RA_BGK_VelocityProfiles.png'.
    %
    % OUTPUT
    %   profiles
    %       Structure containing the reconstructed dimensional profiles,
    %       zero-crossing locations, and flux checks.
    %
    % IMPORTANT
    %   A_o and A_w are liquid-interface velocity amplitudes. Because the
    %   gas is solved kinetically, the macroscopic gas velocity immediately
    %   adjacent to an interface does not generally equal the corresponding
    %   liquid-interface velocity.
    %
    % UNITS
    %   eta                    dimensionless
    %   velocities             m/s
    %   dimensional fluxes     m^2/s
    %   pressure gradient G    Pa/m

    arguments
        result struct
        plotOptions struct = struct()
    end

    %% 1. Default plotting options

    defaults = struct();
    defaults.visible = true;
    defaults.saveFigure = false;
    defaults.outputFile = 'RA_BGK_VelocityProfiles.png';
    defaults.lineWidth = 1.8;

    plotOptions = RA_BGK_MergeStruct(defaults, plotOptions);

    %% 2. Check the model output

    if ~isfield(result, 'success') || ~result.success
        error('RA_BGK:UnsuccessfulResult', ...
              'The supplied model result is not a successful solution.');
    end

    if ~isfield(result, 'solution') || isempty(result.solution)
        error('RA_BGK:MissingSolution', ...
              'The supplied result does not contain a converged solution.');
    end

    state = result.solution;

    requiredFields = { ...
                      'Kn', ...
                      'delta_a', ...
                      'A_o', ...
                      'A_w', ...
                      'G_amplitude'};

    for k = 1:numel(requiredFields)
        if ~isfield(state, requiredFields{k})
            error('RA_BGK:MissingField', ...
                  'The solution does not contain state.%s.', ...
                  requiredFields{k});
        end
    end

    if ~isfield(result, 'options') || ~isfield(result.options, 'kinetic')
        error('RA_BGK:MissingOptions', ...
              'The result does not contain the kinetic solver options.');
    end

    kineticOptions = result.options.kinetic;

    %% 3. Solve the three kinetic basis problems at the converged Kn

    % This calculation reconstructs the full basis profiles at the exact
    % converged Kn. It does not rely on interpolating the stored profiles.
    response = RA_BGK_SolvePlanarResponse( ...
                                          state.Kn, ...
                                          kineticOptions.alpha_oa, ...
                                          kineticOptions.alpha_aw, ...
                                          kineticOptions);

    if ~response.converged
        error('RA_BGK:ProfileSolveFailed', ...
              'The profile reconstruction solve did not converge.');
    end

    if ~isfield(response, 'basis') || numel(response.basis) < 3
        error('RA_BGK:MissingBasisProfiles', ...
              ['RA_BGK_SolvePlanarResponse must return the three basis ', ...
                   'profiles in response.basis.']);
    end

    % Dimensionless normal coordinate.
    nSpace = numel(response.basis{1}.u);
    eta = linspace(0, 1, nSpace).';

    % Dimensionless basis profiles.
    uU0 = response.basis{1}.u;  % unit oil-air interface motion
    uU1 = response.basis{2}.u;  % unit air-water interface motion
    uS  = response.basis{3}.u;  % unit pressure source

    %% 4. Recover the dimensional pressure-source amplitude

    % The interpolation routine stores:
    %
    %       S/G = -2*lambda_mfp^2/(pi*mu_a)
    %
    % Multiplication by G gives a quantity with velocity units.
    if isfield(state, 'kinetic') && ...
            isfield(state.kinetic, 'sourceToPressureGradient')

        sourceToGradient = state.kinetic.sourceToPressureGradient;

    elseif isfield(result, 'fixed') && ...
            isfield(result.fixed, 'lambda_mfp') && ...
            isfield(result.fixed, 'mu_a')

        sourceToGradient = ...
            -2 * result.fixed.lambda_mfp^2 / (pi * result.fixed.mu_a);

    else
        error('RA_BGK:MissingPressureScaling', ...
              'Unable to recover the physical pressure-source scaling.');
    end

    sourceAmplitude = sourceToGradient * state.G_amplitude;

    %% 5. Reconstruct the dimensional velocity components

    % Couette contribution from motion of the two liquid interfaces.
    uOilMotion = state.A_o * uU0;
    uWaterMotion = state.A_w * uU1;

    uCouette = uOilMotion + uWaterMotion;

    % Pressure-driven kinetic contribution.
    uPressure = sourceAmplitude * uS;

    % Total gas velocity in the adopted laboratory/droplet frame.
    uGas = uCouette + uPressure;

    % Velocity relative to the oil-air interface.
    uRelative = uGas - state.A_o;

    % Couette contribution in the oil-interface reference frame.
    uCouetteRelative = uCouette - state.A_o;

    % The pressure component is unchanged by a reference-frame translation.
    uPressureRelative = uPressure;

    %% 6. Calculate dimensional gas fluxes

    h = state.delta_a;

    % Depth-integrated velocity gives flux per unit circumferential width.
    qCouetteRelative = h * trapz(eta, uCouetteRelative);
    qPressure = h * trapz(eta, uPressureRelative);
    qRelative = h * trapz(eta, uRelative);

    % Compare profile integration against the model output.
    if isfield(state, 'q_relative_check')
        qRelativeModel = state.q_relative_check;
    else
        qRelativeModel = NaN;
    end

    %% 7. Locate velocity-reversal points

    etaZeroAbsolute = findZeroCrossings(eta, uGas);
    etaZeroRelative = findZeroCrossings(eta, uRelative);

    % Convert dimensionless eta locations into physical distance from the
    % oil-air interface.
    yZeroAbsolute = h * etaZeroAbsolute;
    yZeroRelative = h * etaZeroRelative;

    %% 8. Calculate cumulative relative flux

    % This shows how positive and negative regions cancel across the layer.
    cumulativeRelativeFlux = h * cumtrapz(eta, uRelative);

    %% 9. Collect outputs

    profiles = struct();

    profiles.eta = eta;
    profiles.y_m = h * eta;
    profiles.y_nm = 1e9 * h * eta;

    profiles.uU0_dimensionless = uU0;
    profiles.uU1_dimensionless = uU1;
    profiles.uS_dimensionless = uS;

    profiles.uOilMotion_m_s = uOilMotion;
    profiles.uWaterMotion_m_s = uWaterMotion;
    profiles.uCouette_m_s = uCouette;
    profiles.uPressure_m_s = uPressure;
    profiles.uGas_m_s = uGas;

    profiles.uCouetteRelative_m_s = uCouetteRelative;
    profiles.uPressureRelative_m_s = uPressureRelative;
    profiles.uRelative_m_s = uRelative;

    profiles.A_o_m_s = state.A_o;
    profiles.A_w_m_s = state.A_w;
    profiles.G_Pa_m = state.G_amplitude;
    profiles.sourceAmplitude_m_s = sourceAmplitude;

    profiles.qCouetteRelative_m2_s = qCouetteRelative;
    profiles.qPressure_m2_s = qPressure;
    profiles.qRelativeIntegrated_m2_s = qRelative;
    profiles.qRelativeModel_m2_s = qRelativeModel;

    profiles.etaZeroAbsolute = etaZeroAbsolute;
    profiles.etaZeroRelative = etaZeroRelative;
    profiles.yZeroAbsolute_nm = 1e9 * yZeroAbsolute;
    profiles.yZeroRelative_nm = 1e9 * yZeroRelative;

    profiles.cumulativeRelativeFlux_m2_s = cumulativeRelativeFlux;

    %% 10. Print numerical summary

    fprintf('\n');
    fprintf('L-BGK GAS-VELOCITY PROFILE RECONSTRUCTION\n');
    fprintf('==========================================\n');
    fprintf('Closure                : %s\n', ...
            string(result.model.gasClosure));
    fprintf('delta_a                 : %.6f nm\n', 1e9 * h);
    fprintf('Kn                      : %.6f\n', state.Kn);
    fprintf('A_o                     : %.6f mm/s\n', 1e3 * state.A_o);
    fprintf('A_w                     : %.6f mm/s\n', 1e3 * state.A_w);
    fprintf('G                       : %.6e Pa/m\n', state.G_amplitude);
    fprintf('Pressure-source amplitude: %.6e m/s\n', sourceAmplitude);
    fprintf('\n');

    fprintf('Relative Couette flux   : %.6e m^2/s\n', ...
            qCouetteRelative);
    fprintf('Pressure-driven flux    : %.6e m^2/s\n', ...
            qPressure);
    fprintf('Integrated net rel. flux: %.6e m^2/s\n', ...
            qRelative);
    fprintf('Model net rel. flux     : %.6e m^2/s\n', ...
            qRelativeModel);
    fprintf('\n');

    if isempty(etaZeroAbsolute)
        fprintf('Absolute gas velocity has no detected zero crossing.\n');
    else
        fprintf('Absolute-velocity zero crossing(s), eta:\n');
        disp(etaZeroAbsolute.');
        fprintf('Absolute-velocity zero crossing(s), y \n');
        disp((1e9 * yZeroAbsolute).');
    end

    if isempty(etaZeroRelative)
        fprintf('Relative gas velocity has no detected zero crossing.\n');
    else
        fprintf('Relative-velocity zero crossing(s), eta:\n');
        disp(etaZeroRelative.');
        fprintf('Relative-velocity zero crossing(s), y \n');
        disp((1e9 * yZeroRelative).');
    end

    %% 11. Create figure

    visibility = 'off';
    if plotOptions.visible
        visibility = 'on';
    end

    figureHandle = figure( ...
                          'Color', 'w', ...
                          'Visible', visibility, ...
                          'Position', [80 80 1250 850]);

    tiledlayout(2, 2, ...
                'Padding', 'compact', ...
                'TileSpacing', 'compact');

    % ----------------------------------------------------------------------
    % Panel A: Absolute gas velocity
    % ----------------------------------------------------------------------

    nexttile;

    plot(eta, 1e3 * uGas, ...
         'k-', ...
         'LineWidth', plotOptions.lineWidth);
    hold on;

    plot(eta, 1e3 * uCouette, ...
         '--', ...
         'LineWidth', 1.4);

    plot(eta, 1e3 * uPressure, ...
         '-.', ...
         'LineWidth', 1.4);

    yline(0, ':', 'LineWidth', 1.0);

    % Plot the liquid-interface velocity amplitudes separately. These points
    % do not imply no slip; they show the velocities of the liquid boundaries.
    plot(0, 1e3 * state.A_o, ...
         'o', ...
         'MarkerSize', 7, ...
         'MarkerFaceColor', [0.85 0.33 0.10], ...
         'MarkerEdgeColor', 'k');

    plot(1, 1e3 * state.A_w, ...
         's', ...
         'MarkerSize', 7, ...
         'MarkerFaceColor', [0 0.45 0.74], ...
         'MarkerEdgeColor', 'k');

    xlabel('\eta=y/\delta_a');
    ylabel('Velocity [mm s^{-1}]');
    title('Gas velocity in the adopted frame');

    legend( ...
           'Total gas velocity', ...
           'Interface-motion contribution', ...
           'Pressure-driven contribution', ...
           'Zero velocity', ...
           'Oil-air interface velocity, A_o', ...
           'Air-water interface velocity, A_w', ...
           'Location', 'best');

    grid on;
    box on;

    % ----------------------------------------------------------------------
    % Panel B: Oil-interface-relative velocity
    % ----------------------------------------------------------------------

    nexttile;

    plot(eta, 1e3 * uRelative, ...
         'k-', ...
         'LineWidth', plotOptions.lineWidth);
    hold on;

    plot(eta, 1e3 * uCouetteRelative, ...
         '--', ...
         'LineWidth', 1.4);

    plot(eta, 1e3 * uPressureRelative, ...
         '-.', ...
         'LineWidth', 1.4);

    yline(0, ':', 'LineWidth', 1.0);

    xlabel('\eta=y/\delta_a');
    ylabel('u_a-A_o [mm s^{-1}]');
    title('Gas velocity relative to the oil-air interface');

    legend( ...
           'Total relative velocity', ...
           'Relative interface-motion contribution', ...
           'Pressure-driven contribution', ...
           'Zero relative velocity', ...
           'Location', 'best');

    grid on;
    box on;

    % ----------------------------------------------------------------------
    % Panel C: Physical distance across the gas film
    % ----------------------------------------------------------------------

    nexttile;

    plot(profiles.y_nm, 1e3 * uRelative, ...
         'k-', ...
         'LineWidth', plotOptions.lineWidth);
    hold on;

    yline(0, ':', 'LineWidth', 1.0);

    if ~isempty(yZeroRelative)
        for k = 1:numel(yZeroRelative)
            xline(1e9 * yZeroRelative(k), ...
                  '--', ...
                  sprintf('Reversal: %.1f nm', 1e9 * yZeroRelative(k)), ...
                  'LabelVerticalAlignment', 'middle');
        end
    end

    xlabel('Distance from oil-air interface [nm]');
    ylabel('u_a-A_o [mm s^{-1}]');
    title('Cross-film relative-flow reversal');

    grid on;
    box on;

    % ----------------------------------------------------------------------
    % Panel D: Cumulative relative flux
    % ----------------------------------------------------------------------

    nexttile;

    plot(eta, cumulativeRelativeFlux, ...
         'Color', [0.49 0.18 0.56], ...
         'LineWidth', plotOptions.lineWidth);
    hold on;

    yline(0, ':', 'LineWidth', 1.0);

    xlabel('\eta=y/\delta_a');
    ylabel('Cumulative relative flux [m^2 s^{-1}]');
    title('Cumulative gas transport across the film');

    grid on;
    box on;

    sgtitle(sprintf( ...
                    'L-BGK velocity reconstruction: %s', ...
                    string(result.model.gasClosure)));

    %% 12. Save the figure if requested

    if plotOptions.saveFigure
        exportgraphics( ...
                       figureHandle, ...
                       plotOptions.outputFile, ...
                       'Resolution', 300);
    end

end

% ==========================================================================
% Local helper: zero-crossing interpolation
% ==========================================================================

function etaZero = findZeroCrossings(eta, velocity)
    % FINDZEROCROSSINGS Locate zero crossings using linear interpolation.

    etaZero = [];

    for k = 1:numel(eta) - 1

        v1 = velocity(k);
        v2 = velocity(k + 1);

        if ~isfinite(v1) || ~isfinite(v2)
            continue
        end

        if v1 == 0
            etaZero(end + 1, 1) = eta(k);

        elseif v1 * v2 < 0

            fraction = -v1 / (v2 - v1);

            etaCross = eta(k) + ...
                fraction * (eta(k + 1) - eta(k));

            etaZero(end + 1, 1) = etaCross;
        end
    end

    if velocity(end) == 0
        etaZero(end + 1, 1) = eta(end);
    end

    % Remove repeated zero crossings produced by exact grid-node zeros.
    if ~isempty(etaZero)
        etaZero = uniquetol(etaZero, 1e-10);
    end

end
