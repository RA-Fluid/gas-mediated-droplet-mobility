function output = RA_BGK_PlotSphericalVelocityReconstruction( ...
                                                             profilesClosed, profilesOpen, plotOptions)
    % RA_BGK_PLOTSPHERICALVELOCITYRECONSTRUCTION
    % Reconstruct and plot the tangential L-BGK gas velocity on a half sphere.
    %
    % PURPOSE
    %   Maps the locally planar, cross-film gas-velocity profiles onto the
    %   leading spherical harmonic:
    %
    %       u_theta(eta,theta) = u_relative(eta)*sin(theta)
    %
    %   where:
    %
    %       eta   = (r-R_0)/delta_a
    %       theta = polar angle, 0 <= theta <= pi
    %
    %   The function compares:
    %
    %       1. zero-relative-flux closure;
    %       2. zero-pressure-gradient closure.
    %
    %   Colours show the signed tangential velocity relative to the oil-air
    %   interface. Arrows show the local tangential direction and relative
    %   magnitude.
    %
    % IMPORTANT PHYSICAL LIMITATION
    %   This figure is a reconstructed tangential velocity field based on:
    %
    %       - a uniform gas-layer thickness;
    %       - the locally planar L-BGK profile;
    %       - the leading sin(theta) spherical harmonic.
    %
    %   The calculation does not resolve radial gas velocity, global
    %   streamlines, or a two-dimensional recirculation cell. The output
    %   should therefore not be described as a vortex or streamline solution.
    %
    % REQUIRED INPUT STRUCTURES
    %   profilesClosed and profilesOpen are outputs from:
    %
    %       RA_BGK_PlotVelocityProfiles
    %
    %   Each structure must contain:
    %
    %       .eta
    %       .uRelative_m_s
    %       .uGas_m_s
    %       .A_o_m_s
    %       .A_w_m_s
    %       .G_Pa_m
    %       .qRelativeIntegrated_m2_s
    %       .etaZeroRelative
    %
    % OPTIONAL SETTINGS
    %   plotOptions.displayThickness
    %       Displayed shell thickness relative to the inner radius.
    %       Default: 0.18.
    %
    %   plotOptions.nTheta
    %       Number of polar-angle points. Default: 181.
    %
    %   plotOptions.nEta
    %       Number of cross-film plotting points. Default: 81.
    %
    %   plotOptions.nThetaArrows
    %       Number of angular arrow positions. Default: 15.
    %
    %   plotOptions.nEtaArrows
    %       Number of radial arrow positions. Default: 7.
    %
    %   plotOptions.arrowScale
    %       Visual arrow-length multiplier. Default: 0.11.
    %
    %   plotOptions.showAbsoluteVelocity
    %       false plots velocity relative to the oil-air interface.
    %       true plots gas velocity in the adopted absolute frame.
    %       Default: false.
    %
    %   plotOptions.saveFigure
    %       Save the generated figure. Default: false.
    %
    %   plotOptions.outputFile
    %       Output image filename.
    %
    % OUTPUT
    %   output
    %       Structure containing plotting grids, mapped velocity fields,
    %       reversal locations, and figure handle.
    %
    % EXAMPLE
    %   figureOutput = RA_BGK_PlotSphericalVelocityReconstruction( ...
    %       profilesClosed,profilesOpen);
    %
    %   The shell thickness is exaggerated only geometrically. The displayed
    %   velocity values and directions retain the calculated L-BGK results.

    arguments
        profilesClosed struct
        profilesOpen struct
        plotOptions struct = struct()
    end

    %% 1. Default plotting options

    defaults = struct();

    % The physical gas layer is too thin to see at the true geometric scale.
    % This parameter controls only the displayed annular thickness.
    defaults.displayThickness = 0.18;

    defaults.nTheta = 181;
    defaults.nEta = 81;

    defaults.nThetaArrows = 15;
    defaults.nEtaArrows = 7;

    % Controls arrow length in displayed-radius units.
    defaults.arrowScale = 0.11;

    defaults.showAbsoluteVelocity = false;

    defaults.visible = true;
    defaults.saveFigure = false;
    defaults.outputFile = ...
        'RA_BGK_SphericalVelocityReconstruction.png';

    defaults.colormapName = 'turbo';
    defaults.contourLevels = 101;
    defaults.showReversalSurface = true;
    defaults.showInterfaceLabels = true;

    plotOptions = RA_BGK_MergeStruct(defaults, plotOptions);

    %% 2. Validate the supplied profile structures

    validateProfile(profilesClosed, 'profilesClosed');
    validateProfile(profilesOpen, 'profilesOpen');

    if plotOptions.displayThickness <= 0
        error('displayThickness must be positive.');
    end

    if plotOptions.nTheta < 21
        error('nTheta must be at least 21.');
    end

    if plotOptions.nEta < 11
        error('nEta must be at least 11.');
    end

    %% 3. Select the velocity field to display

    if plotOptions.showAbsoluteVelocity

        velocityClosed = profilesClosed.uGas_m_s;
        velocityOpen = profilesOpen.uGas_m_s;

        velocityLabel = 'Tangential gas velocity [mm s^{-1}]';
        frameDescription = 'Gas velocity in the adopted frame';

    else

        velocityClosed = profilesClosed.uRelative_m_s;
        velocityOpen = profilesOpen.uRelative_m_s;

        velocityLabel = ...
            'Tangential velocity relative to oil-air interface [mm s^{-1}]';

        frameDescription = ...
            'Gas velocity relative to the oil-air interface';

    end

    %% 4. Build a common cross-film and angular grid

    etaPlot = linspace(0, 1, plotOptions.nEta).';

    thetaPlot = linspace(0, pi, plotOptions.nTheta);

    % Interpolate the original L-BGK profiles onto the plotting grid.
    uClosed = interp1( ...
                      profilesClosed.eta(:), ...
                      velocityClosed(:), ...
                      etaPlot, ...
                      'pchip');

    uOpen = interp1( ...
                    profilesOpen.eta(:), ...
                    velocityOpen(:), ...
                    etaPlot, ...
                    'pchip');

    % Convert to millimetres per second for plotting.
    uClosed_mm_s = 1e3 * uClosed;
    uOpen_mm_s = 1e3 * uOpen;

    %% 5. Construct the exaggerated displayed shell geometry

    % The inner displayed radius is normalised to unity.
    displayInnerRadius = 1;

    % The displayed gas-layer thickness is exaggerated.
    displayOuterRadius = ...
        displayInnerRadius + plotOptions.displayThickness;

    radiusPlot = ...
        displayInnerRadius + ...
        plotOptions.displayThickness * etaPlot;

    % Mesh arrays have size nEta x nTheta.
    [ThetaGrid, RadiusGrid] = meshgrid(thetaPlot, radiusPlot);

    % Meridional-plane coordinates:
    %
    %       x = r*sin(theta)
    %       z = r*cos(theta)
    %
    % The half sphere is therefore displayed in the x >= 0 plane.
    X = RadiusGrid .* sin(ThetaGrid);
    Z = RadiusGrid .* cos(ThetaGrid);

    %% 6. Map the cross-film profiles onto the spherical harmonic

    % The leading spherical angular dependence is sin(theta).
    %
    % u_theta(eta,theta) = u_profile(eta)*sin(theta)
    %
    % Each column corresponds to one theta value.
    UthetaClosed = uClosed_mm_s .* sin(thetaPlot);
    UthetaOpen = uOpen_mm_s .* sin(thetaPlot);

    %% 7. Convert tangential velocity into Cartesian vector components

    % Unit vector in the increasing-theta direction:
    %
    %       e_theta = cos(theta)*e_x - sin(theta)*e_z
    %
    % Therefore:
    %
    %       U_x = U_theta*cos(theta)
    %       U_z = -U_theta*sin(theta)

    UxClosed = UthetaClosed .* cos(thetaPlot);
    UzClosed = -UthetaClosed .* sin(thetaPlot);

    UxOpen = UthetaOpen .* cos(thetaPlot);
    UzOpen = -UthetaOpen .* sin(thetaPlot);

    %% 8. Use one symmetric colour range for both closures

    maxAbsVelocity = max(abs([UthetaClosed(:); UthetaOpen(:)]));

    if maxAbsVelocity <= 0 || ~isfinite(maxAbsVelocity)
        maxAbsVelocity = 1;
    end

    colourLimits = [-maxAbsVelocity, maxAbsVelocity];

    %% 9. Construct a sparse vector grid

    etaArrow = linspace(0.08, 0.92, plotOptions.nEtaArrows).';
    thetaArrow = linspace(0.08 * pi, 0.92 * pi, ...
                          plotOptions.nThetaArrows);

    radiusArrow = ...
        displayInnerRadius + ...
        plotOptions.displayThickness * etaArrow;

    [ThetaArrow, RadiusArrow] = meshgrid(thetaArrow, radiusArrow);

    XArrow = RadiusArrow .* sin(ThetaArrow);
    ZArrow = RadiusArrow .* cos(ThetaArrow);

    uClosedArrow = interp1( ...
                           etaPlot, uClosed_mm_s, etaArrow, 'pchip');

    uOpenArrow = interp1( ...
                         etaPlot, uOpen_mm_s, etaArrow, 'pchip');

    UthetaClosedArrow = uClosedArrow .* sin(thetaArrow);
    UthetaOpenArrow = uOpenArrow .* sin(thetaArrow);

    UxClosedArrow = UthetaClosedArrow .* cos(thetaArrow);
    UzClosedArrow = -UthetaClosedArrow .* sin(thetaArrow);

    UxOpenArrow = UthetaOpenArrow .* cos(thetaArrow);
    UzOpenArrow = -UthetaOpenArrow .* sin(thetaArrow);

    % Normalise vector lengths for readability while preserving direction.
    [UxClosedDisplay, UzClosedDisplay] = normaliseArrows( ...
                                                        UxClosedArrow, UzClosedArrow, ...
                                                        maxAbsVelocity, plotOptions.arrowScale * 3);

    [UxOpenDisplay, UzOpenDisplay] = normaliseArrows( ...
                                                    UxOpenArrow, UzOpenArrow, ...
                                                    maxAbsVelocity, plotOptions.arrowScale);

    %% 10. Determine the cross-film reversal locations

    etaZeroClosed = findZeroCrossings(etaPlot, uClosed_mm_s);
    etaZeroOpen = findZeroCrossings(etaPlot, uOpen_mm_s);

    %% 11. Create the figure

    visibility = 'off';

    if plotOptions.visible
        visibility = 'on';
    end

    figureHandle = figure( ...
                          'Color', 'w', ...
                          'Visible', visibility, ...
                          'Position', [70 90 1450 690]);

    layout = tiledlayout(1, 2, ...
                         'Padding', 'compact', ...
                         'TileSpacing', 'compact');

    % ----------------------------------------------------------------------
    % Panel A: zero-relative-flux closure
    % ----------------------------------------------------------------------

    axisClosed = nexttile;

    plotSphericalPanel( ...
                       axisClosed, ...
                       X, Z, UthetaClosed, ...
                       XArrow, ZArrow, UxClosedDisplay, UzClosedDisplay, ...
                       etaZeroClosed, ...
                       displayInnerRadius, displayOuterRadius, ...
                       plotOptions, ...
                       colourLimits, ...
                       'Zero-relative-flux closure', ...
                       sprintf(['G = %.2e Pa m^{-1}, ', ...
                 'q_{rel} = %.2e m^2 s^{-1}'], ...
                profilesClosed.G_Pa_m, ...
                profilesClosed.qRelativeIntegrated_m2_s));

    % ----------------------------------------------------------------------
    % Panel B: zero-pressure-gradient closure
    % ----------------------------------------------------------------------

    axisOpen = nexttile;

    plotSphericalPanel( ...
                       axisOpen, ...
                       X, Z, UthetaOpen, ...
                       XArrow, ZArrow, UxOpenDisplay, UzOpenDisplay, ...
                       etaZeroOpen, ...
                       displayInnerRadius, displayOuterRadius, ...
                       plotOptions, ...
                       colourLimits, ...
                       'Pressure-gradient-free closure', ...
                       sprintf(['G = %.2e Pa m^{-1}, ', ...
                 'q_{rel} = %.2e m^2 s^{-1}'], ...
                profilesOpen.G_Pa_m, ...
                profilesOpen.qRelativeIntegrated_m2_s));

    %% 12. Apply one colormap and add a shared colour bar

    % A tiledlayout object is not a valid first argument to colormap in some
    % MATLAB releases. Apply the colormap to the parent figure instead.
    %
    % Convert the requested colormap name into an explicit RGB matrix. This
    % also avoids version-dependent handling of string-valued colormap names.
    try
        colourMapValues = feval( ...
                                plotOptions.colormapName, ...
                                256);
    catch ME
        error('RA_BGK:InvalidColormap', ...
              'Unable to construct the requested colormap "%s": %s', ...
              string(plotOptions.colormapName), ...
              ME.message);
    end

    colormap(figureHandle, colourMapValues);

    % Apply identical signed velocity limits to both panels, so that the same
    % colour represents the same velocity under both gas-transport closures.
    clim(axisClosed, colourLimits);
    clim(axisOpen, colourLimits);

    % Create the colour bar from one of the axes. Because both axes use the
    % same figure colormap and identical colour limits, this colour bar applies
    % to both panels.
    colourBar = colorbar(axisOpen);

    % Newer MATLAB releases allow a colour bar to occupy the east tile of a
    % tiled layout. Retain a fallback for releases that do not support this
    % property.
    try
        colourBar.Layout.Tile = 'east';
    catch
        % Fallback placement for older MATLAB releases.
        drawnow;

        colourBar.Position = [ ...
                              0.925, ... % left
                              0.18,  ... % bottom
                              0.018, ... % width
                              0.67];     % height
    end

    colourBar.Label.String = velocityLabel;
    colourBar.Label.Interpreter = 'tex';

    % Add the overall figure title.
    sgtitle(layout, ...
            sprintf([ ...
                 'Reconstructed tangential gas velocity on a uniform ', ...
                 'spherical shell\n%s; shell thickness exaggerated'], ...
                frameDescription), ...
            'FontWeight', 'bold');
    %% 13. Save the figure if requested

    if plotOptions.saveFigure

        exportgraphics( ...
                       figureHandle, ...
                       plotOptions.outputFile, ...
                       'Resolution', 300);

    end

    %% 14. Return the mapped fields

    output = struct();

    output.figureHandle = figureHandle;
    output.axisClosed = axisClosed;
    output.axisOpen = axisOpen;

    output.eta = etaPlot;
    output.theta = thetaPlot;
    output.displayRadius = radiusPlot;

    output.X = X;
    output.Z = Z;

    output.uClosedCrossFilm_mm_s = uClosed_mm_s;
    output.uOpenCrossFilm_mm_s = uOpen_mm_s;

    output.UthetaClosed_mm_s = UthetaClosed;
    output.UthetaOpen_mm_s = UthetaOpen;

    output.UxClosed_mm_s = UxClosed;
    output.UzClosed_mm_s = UzClosed;

    output.UxOpen_mm_s = UxOpen;
    output.UzOpen_mm_s = UzOpen;

    output.etaZeroClosed = etaZeroClosed;
    output.etaZeroOpen = etaZeroOpen;

    output.colourLimits_mm_s = colourLimits;
    output.displayThickness = plotOptions.displayThickness;

end

% ==========================================================================
% Local function: draw one half-spherical velocity panel
% ==========================================================================

function plotSphericalPanel( ...
                            axisHandle, ...
                            X, Z, Utheta, ...
                            XArrow, ZArrow, UxArrow, UzArrow, ...
                            etaZero, ...
                            innerRadius, outerRadius, ...
                            options, colourLimits, ...
                            panelTitle, panelSubtitle)

    axes(axisHandle); %#ok<LAXES>

    % Signed tangential velocity contour.
    surfaceHandle = surf( ...
                         axisHandle, ...
                         X, Z, zeros(size(X)), Utheta, ...
                         'EdgeColor', 'none', ...
                         'FaceColor', 'interp');

    view(axisHandle, 2);
    hold(axisHandle, 'on');

    % Tangential vectors.
    quiver( ...
           axisHandle, ...
           XArrow, ZArrow, ...
           UxArrow, UzArrow, ...
           0, ...
           'k', ...
           'LineWidth', 0.9, ...
           'MaxHeadSize', 0.7);

    % Inner and outer shell boundaries.
    thetaBoundary = linspace(0, pi, 400);

    plot( ...
         axisHandle, ...
         innerRadius * sin(thetaBoundary), ...
         innerRadius * cos(thetaBoundary), ...
         'k-', ...
         'LineWidth', 1.6);

    plot( ...
         axisHandle, ...
         outerRadius * sin(thetaBoundary), ...
         outerRadius * cos(thetaBoundary), ...
         'k-', ...
         'LineWidth', 1.6);

    % Plot reversal surfaces.
    if options.showReversalSurface && ~isempty(etaZero)

        for k = 1:numel(etaZero)

            reversalRadius = ...
                innerRadius + ...
                options.displayThickness * etaZero(k);

            plot( ...
                 axisHandle, ...
                 reversalRadius * sin(thetaBoundary), ...
                 reversalRadius * cos(thetaBoundary), ...
                 'w--', ...
                 'LineWidth', 2.3);

            plot( ...
                 axisHandle, ...
                 reversalRadius * sin(thetaBoundary), ...
                 reversalRadius * cos(thetaBoundary), ...
                 'k--', ...
                 'LineWidth', 0.9);
        end
    end

    % Centreline separating the visualised half-plane.
    plot( ...
         axisHandle, ...
         [0 0], ...
         [-outerRadius outerRadius], ...
         'Color', [0.35 0.35 0.35], ...
         'LineStyle', ':', ...
         'LineWidth', 0.8);

    axis(axisHandle, 'equal');
    axis(axisHandle, 'off');

    xlim(axisHandle, [0, 1.12 * outerRadius]);
    ylim(axisHandle, [-1.08 * outerRadius, 1.08 * outerRadius]);

    clim(axisHandle, colourLimits);

    title( ...
          axisHandle, ...
          {panelTitle, panelSubtitle}, ...
          'FontWeight', 'bold');

    % Add labels only when requested.
    if options.showInterfaceLabels

        text( ...
             axisHandle, ...
             0.08 * outerRadius, ...
             0.89 * innerRadius, ...
             'Oil-air interface', ...
             'FontSize', 9, ...
             'Rotation', -75, ...
             'HorizontalAlignment', 'left');

        text( ...
             axisHandle, ...
             0.13 * outerRadius, ...
             1.02 * outerRadius, ...
             'Air-water interface', ...
             'FontSize', 9, ...
             'Rotation', -75, ...
             'HorizontalAlignment', 'left');

        if options.showReversalSurface && ~isempty(etaZero)

            labelRadius = ...
                innerRadius + ...
                options.displayThickness * etaZero(1);

            text( ...
                 axisHandle, ...
                 labelRadius, ...
                 0, ...
                 sprintf('\\eta^* = %.3f', etaZero(1)), ...
                 'Color', 'k', ...
                 'BackgroundColor', 'w', ...
                 'Margin', 2, ...
                 'FontSize', 9, ...
                 'HorizontalAlignment', 'center');
        end
    end

    % Maintain handle so MATLAB does not optimise away the surface.
    if isempty(surfaceHandle)
        error('Unable to create the spherical velocity surface.');
    end

end

% ==========================================================================
% Local function: normalise arrows for visual display
% ==========================================================================

function [UxDisplay, UzDisplay] = normaliseArrows( ...
                                                 Ux, Uz, maxVelocity, arrowScale)

    magnitude = sqrt(Ux.^2 + Uz.^2);

    % Preserve relative velocity magnitude, but express arrow lengths in the
    % normalised geometric coordinates used by the displayed sphere.
    scaleFactor = arrowScale / max(maxVelocity, realmin);

    UxDisplay = scaleFactor * Ux;
    UzDisplay = scaleFactor * Uz;

    % Suppress numerical noise in nearly zero vectors.
    small = magnitude < 1e-12 * max(maxVelocity, realmin);

    UxDisplay(small) = 0;
    UzDisplay(small) = 0;

end

% ==========================================================================
% Local function: validate a reconstructed profile structure
% ==========================================================================

function validateProfile(profile, inputName)

    requiredFields = { ...
                      'eta', ...
                      'uRelative_m_s', ...
                      'uGas_m_s', ...
                      'A_o_m_s', ...
                      'A_w_m_s', ...
                      'G_Pa_m', ...
                      'qRelativeIntegrated_m2_s'};

    for k = 1:numel(requiredFields)

        if ~isfield(profile, requiredFields{k})

            error('RA_BGK:MissingProfileField', ...
                  '%s does not contain the field %s.', ...
                  inputName, requiredFields{k});
        end
    end

    if numel(profile.eta) ~= numel(profile.uRelative_m_s)

        error('RA_BGK:ProfileSizeMismatch', ...
              '%s.eta and %s.uRelative_m_s have different lengths.', ...
              inputName, inputName);
    end

    if numel(profile.eta) ~= numel(profile.uGas_m_s)

        error('RA_BGK:ProfileSizeMismatch', ...
              '%s.eta and %s.uGas_m_s have different lengths.', ...
              inputName, inputName);
    end

end

% ==========================================================================
% Local function: locate cross-film zero crossings
% ==========================================================================

function etaZero = findZeroCrossings(eta, velocity)

    etaZero = [];

    for k = 1:numel(eta) - 1

        v1 = velocity(k);
        v2 = velocity(k + 1);

        if ~isfinite(v1) || ~isfinite(v2)
            continue
        end

        if abs(v1) <= 10 * eps(max(1, max(abs(velocity))))

            etaZero(end + 1, 1) = eta(k); %#ok<AGROW>

        elseif v1 * v2 < 0

            fraction = -v1 / (v2 - v1);

            etaCross = ...
                eta(k) + ...
                fraction * (eta(k + 1) - eta(k));

            etaZero(end + 1, 1) = etaCross; %#ok<AGROW>
        end
    end

    if abs(velocity(end)) <= 10 * eps(max(1, max(abs(velocity))))
        etaZero(end + 1, 1) = eta(end);
    end

    if ~isempty(etaZero)
        etaZero = uniquetol(etaZero, 1e-10);
    end

end
