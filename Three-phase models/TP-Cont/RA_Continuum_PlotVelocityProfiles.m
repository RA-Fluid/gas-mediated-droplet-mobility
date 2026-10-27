function profiles = RA_Continuum_PlotVelocityProfiles(resultClosed, resultOpen, plotOptions)
    % RA_CONTINUUM_PLOTVELOCITYPROFILES Compare the two continuum profiles.
    arguments
        resultClosed struct
        resultOpen struct
        plotOptions struct = struct()
    end
    def = struct('visible', true, 'saveFigure', false, ...
               'outputFile', 'RA_Continuum_ClosureProfiles.png', 'nPoints', 201);
    plotOptions = RA_Continuum_MergeStruct(def, plotOptions);
    c = RA_Continuum_ReconstructVelocityProfile(resultClosed, plotOptions.nPoints);
    o = RA_Continuum_ReconstructVelocityProfile(resultOpen, plotOptions.nPoints);
    vis = 'off';
    if plotOptions.visible
        vis = 'on';
    end
    fig = figure('Color', 'w', 'Visible', vis, 'Position', [80 100 1250 500]);
    tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile;
    plot(c.eta, 1e3 * c.uRelative_m_s, 'k-', 'LineWidth', 1.8);
    hold on;
    plot(c.eta, 1e3 * c.uCouetteRelative_m_s, '--', 'LineWidth', 1.4);
    plot(c.eta, 1e3 * c.uPressure_m_s, '-.', 'LineWidth', 1.4);
    yline(0, ':');
    xlabel('eta');
    ylabel('u_a-A_o [mm/s]');
    title('Zero relative flux');
    legend('Total', 'Interface motion', 'Pressure response', 'Zero', 'Location', 'best');
    grid on;
    nexttile;
    plot(o.eta, 1e3 * o.uRelative_m_s, 'k-', 'LineWidth', 1.8);
    hold on;
    plot(o.eta, 1e3 * o.uCouetteRelative_m_s, '--', 'LineWidth', 1.4);
    plot(o.eta, 1e3 * o.uPressure_m_s, '-.', 'LineWidth', 1.4);
    yline(0, ':');
    xlabel('eta');
    ylabel('u_a-A_o [mm/s]');
    title('Zero pressure gradient');
    legend('Total', 'Interface motion', 'Pressure response', 'Zero', 'Location', 'best');
    grid on;
    if plotOptions.saveFigure
        exportgraphics(fig, plotOptions.outputFile, 'Resolution', 300);
    end
    profiles = struct('closed', c, 'open', o, 'figureHandle', fig);
end
