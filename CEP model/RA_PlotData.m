%% initialise

clc
close all
clear

% Set default font to Helvetica for all new figures
set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultLegendFontName', 'Arial');

%% Some main plot settings

% oil name
oil_nam = 'SO-20cst' ;

% plot magnifier coefficient
plotMag = 3 ;
% base figure width
baseWdith = 8.9 ; % [mm]
% Aspect Ratio (Width / Height)
aspect_ratio = 4/3 ;

%% saturation levels

Si_W = (1:0.05:1.3)' ;
Si_O = (1:0.05:1.3)' ;

nW = length(Si_W) ;
nO = length(Si_O) ;

% datafile names
datFilNam = cell(nW, nO) ;
for i = 1 : nW
    for j = 1 : nO

        baseName = sprintf('RA_Output_Si_W_%.2f_Si_O_%.2f.mat', Si_W(i), Si_O(j)) ;
        datFilNam{i,j} = fullfile('output', baseName);

    end
end

% selected files for plotting
filNam_PlotA = datFilNam(:, 1) ;
filNam_PlotB = datFilNam(1,:) ;

%% Read the data for Hexadecane

%----------------------------------------------------------------------
% Oil sat. constant, Water sat. changes
%----------------------------------------------------------------------
% selected diameter to plot the changes of air layer thickness
Dd_select = 250 ; % [um]
% the oil chooser id
idOil = 1 ; % hexadecane

delta_nm_Hex_SwVaries = cell(nW, 1) ;
delta_Vis_Hex_SwVaries = zeros() ;
for i = 1 : nW

    % read the full data for one saturation condition
    dat = load(filNam_PlotA{i}) ;

    % number of droplets
    nD = size(dat.RA_Output, 1) ;
    % loop over the droplet sizes
    for j = 1 : nD

        % store the calcualted thickness at the visualisation height
        delta_Vis_Hex_SwVaries(j, i) = dat.RA_Output(j, idOil).delta_inf ; % [nm]

        % read the droplet size
        Dd = dat.RA_Output(j, idOil).datInput.Dd * 1e6 ; % [um]
        % save the temporal thickness data only for the prescribed diamater
        if abs( Dd - Dd_select) < 1
            a = dat.RA_Output(j, idOil).t ;
            b = dat.RA_Output(j, idOil).delta_nm ;
            delta_nm_Hex_SwVaries{i,1} = [a, b] ; % [nm]
        end

    end

end

%----------------------------------------------------------------------
% Oil sat. varies, Water sat. constant
%----------------------------------------------------------------------
% the oil chooser id
idOil = 1 ; % hexadecane

delta_nm_Hex_SoVaries = cell(nO, 1) ;
delta_Vis_Hex_SoVaries = zeros() ;
for i = 1 : nO

    % read the full data for one saturation condition
    dat = load(filNam_PlotB{i}) ;

    % number of droplets
    nD = size(dat.RA_Output, 1) ;
    % loop over the droplet sizes
    for j = 1 : nD

        % store the calcualted thickness at the visualisation height
        delta_Vis_Hex_SoVaries(j, i) = dat.RA_Output(j, idOil).delta_inf ; % [nm]

        % read the droplet size
        Dd = dat.RA_Output(j, idOil).datInput.Dd * 1e6 ; % [um]
        % save the temporal thickness data only for the prescribed diamater
        if abs( Dd - Dd_select) < 1
            a = dat.RA_Output(j, idOil).t ;
            b = dat.RA_Output(j, idOil).delta_nm ;
            delta_nm_Hex_SoVaries{i,1} = [a, b] ; % [nm]
        end

    end

end

%% Read the data for 20 cst PDMS
%----------------------------------------------------------------------
% Oil sat. constant, Water sat. changes
%----------------------------------------------------------------------
% selected diameter to plot the changes of air layer thickness
Dd_select = 250 ; % [um]
% the oil chooser id
idOil = 2 ; % 20 cst PDMS

delta_nm_20cst_SwVaries = cell(nW, 1) ;
delta_Vis_20cst_SwVaries = zeros() ;
for i = 1 : nW

    % read the full data for one saturation condition
    dat = load(filNam_PlotA{i}) ;

    % number of droplets
    nD = size(dat.RA_Output, 1) ;
    % loop over the droplet sizes
    for j = 1 : nD

        % store the calcualted thickness at the visualisation height
        delta_Vis_20cst_SwVaries(j, i) = dat.RA_Output(j, idOil).delta_inf ; % [nm]

        % read the droplet size
        Dd = dat.RA_Output(j, idOil).datInput.Dd * 1e6 ; % [um]
        % save the temporal thickness data only for the prescribed diamater
        if abs( Dd - Dd_select) < 1
            a = dat.RA_Output(j, idOil).t ;
            b = dat.RA_Output(j, idOil).delta_nm ;
            delta_nm_20cst_SwVaries{i,1} = [a, b] ; % [nm]
        end

    end

end

%----------------------------------------------------------------------
% Oil sat. varies, Water sat. constant
%----------------------------------------------------------------------
% the oil chooser id
idOil = 2 ; % 20 cst PDMS

delta_nm_20cst_SoVaries = cell(nO, 1) ;
delta_Vis_20cst_SoVaries = zeros() ;
for i = 1 : nO

    % read the full data for one saturation condition
    dat = load(filNam_PlotB{i}) ;

    % number of droplets
    nD = size(dat.RA_Output, 1) ;
    % loop over the droplet sizes
    for j = 1 : nD

        % store the calcualted thickness at the visualisation height
        delta_Vis_20cst_SoVaries(j, i) = dat.RA_Output(j, idOil).delta_inf ; % [nm]

        % read the droplet size
        Dd = dat.RA_Output(j, idOil).datInput.Dd * 1e6 ; % [um]
        % save the temporal thickness data only for the prescribed diamater
        if abs( Dd - Dd_select) < 1
            a = dat.RA_Output(j, idOil).t ;
            b = dat.RA_Output(j, idOil).delta_nm ;
            delta_nm_20cst_SoVaries{i,1} = [a, b] ; % [nm]
        end
    end

end

%% Figure A: Hexadecane | water saturation varies

is_plt = 0 ;
if isequal( is_plt, 1 )
    %----------------------------------------------------------------------
    % plot properties
    %----------------------------------------------------------------------
    % number of data sets
    nDat = 7 ;
    % dataset collection
    dat_plt.dat = delta_nm_Hex_SwVaries(end:-1:1) ;
    % symbolsOS_PU
    dat_plt.Symbols = { '-', '--', '-.', '-', '--', '-.', '-' } ;
    % line colours
    dat_plt.Colors = {RA_getColorHex('Teal'), RA_getColorHex('DarkCyan'), RA_getColorHex('CadetBlue'), ...
        RA_getColorHex('SteelBlue'), RA_getColorHex('LightSeaGreen'),...
        RA_getColorHex('MediumAquamarine'), RA_getColorHex('DarkSlateGray')} ;
    % linewidths
    c = cell(1, nDat) ;
    c(:) = {0.75*plotMag} ;
    dat_plt.LineWidths = c ;
    % line alpha
    c = cell(1, nDat) ;
    c(:) = {1} ;
    dat_plt.LineAlpha = c ;

    % marker sizes
    c = cell(1, nDat) ;
    c(:) = {5*plotMag} ;
    dat_plt.MarkerSizes = c ;
    % edge colours
    dat_plt.MarkerEdgeColors = dat_plt.Colors ;
    % marker face colours
    c = cell(1, nDat) ;
    c(:) = {'none'} ;
    dat_plt.MarkerFaceColors = c ;
    % marker face alpha
    c = cell(1, nDat) ;
    c(:) = {1};
    dat_plt.MarkerAlpha = c ;

    % aspect ratio of the data
    dat_plt.ASR = 1/aspect_ratio ;
    % width of the figure
    dat_plt.figWidth_cm = baseWdith *plotMag ; % [cm]
    % fontsize
    dat_plt.AxesFontSize = 8 *plotMag ;
    % interpreter
    dat_plt.interpreter = 'tex' ;
    % xlabel
    dat_plt.XLabel = '{\itt}{\fontsize{8} }(s)' ;
    % ylabel
    dat_plt.YLabel = '{\it\delta}_{a}{\fontsize{8} }(nm)' ;
    % x limit range
    a = dat_plt.dat{1} ;
    b = dat_plt.dat{end} ;
    dat_plt.XLim = [1e-4, 1.1*max(a(:,1))] ;
    % y limit range
    dat_plt.YLim = [0, 150] ; %[0.9*min(b(:,2)), 1.1*max(b(:,2))] ;
    % number of x-axis labels [negative, positive]
    dat_plt.XLabNum = [1, 3] ;
    dat_plt.XTickValues = [1e-4, 1e-3, 1e-2, 0.1, 1, 10 ] ; %, max(a(:,1))
    dat_plt.XTickLabels = {'10^{-4}', '10^{-3}', '10^{-2}', '0.1', '1', '10'} ; %, num2str( max(a(:,1)), '%.0f' )
    % number of y-axis labels [negative, positive]
    dat_plt.YLabNum = [1, 5] ;
    % x-aixs labels' format
    dat_plt.XLabForm = '%.0f' ;
    % y-aixs labels' format
    dat_plt.YLabForm = '%.0f' ;
    % x-axis scale type
    dat_plt.XScale = 'log' ;
    % y-axis scale type
    dat_plt.YScale = 'linear' ;
    % axis direction
    dat_plt.AxDir = 'left' ;
    % type of axis
    dat_plt.axType = 'single' ;

    dat_plt.splin = 'no' ;
    dat_plt.er_num = 3 ;

    % plot the data
    RA_Plot_1D(dat_plt) ;
    grid on

    % legend
    leg_tex = cell(nW, 1) ;
    for i = 1 : nW

        f = 1 + 0.05 * (i-1) ;
        leg_tex{i} = num2str(f, '%.2f') ;

    end
    leg_tex = leg_tex(end:-1:1) ;
    legend(leg_tex, 'fontsize', dat_plt.AxesFontSize, 'interpreter', 'tex', 'location', 'northwest', 'NumColumns', 1 ) ;

    % add text
    % text(0.5, 120, ['10^3\langle\beta\rangle{\fontsize{3} }={\fontsize{3} }', num2str(beta*1e3, '%.0f'), '\pm', num2str(sig_beta*1e3, '%.0f') ], 'interpreter', 'tex', 'fontsize', dat_plt.AxesFontSize ) ;
    % text(1.2*min(DRange), 1.1*min(URange), ['{\itD}_d{\fontsize{3} }={\fontsize{3} }', num2str(DRange(1), '%.0f'), '-', num2str(DRange(2), '%.0f'), '{\fontsize{3} }', char(181), 'm' ], 'interpreter', 'tex', 'fontsize', dat_plt.AxesFontSize ) ;

    %----------------------------------------------------------------------
    % Export the figure as a vectorised PDF
    %----------------------------------------------------------------------
    exportgraphics( gca, 'Hexadecane-SiO_1.00_SiW_varies.svg', 'ContentType', 'vector') ;
end
clear dat_plt
close

%% Figure B: 20 cst | water saturation varies

is_plt = 1 ;
if isequal( is_plt, 1 )
    %----------------------------------------------------------------------
    % plot properties
    %----------------------------------------------------------------------
    % number of data sets
    nDat = 7 ;
    % dataset collection
    dat_plt.dat = delta_nm_20cst_SwVaries(end:-1:1) ;
    % symbolsOS_PU
    dat_plt.Symbols = { '-', '--', '-.', '-', '--', '-.', '-' } ;
    % line colours
    dat_plt.Colors = {RA_getColorHex('DarkRed'), RA_getColorHex('FireBrick'), RA_getColorHex('Crimson'), ...
        RA_getColorHex('IndianRed'), RA_getColorHex('Tomato'),...
        RA_getColorHex('OrangeRed'), RA_getColorHex('DarkOrange')} ;
    % linewidths
    c = cell(1, nDat) ;
    c(:) = {0.75*plotMag} ;
    dat_plt.LineWidths = c ;
    % line alpha
    c = cell(1, nDat) ;
    c(:) = {1} ;
    dat_plt.LineAlpha = c ;

    % marker sizes
    c = cell(1, nDat) ;
    c(:) = {5*plotMag} ;
    dat_plt.MarkerSizes = c ;
    % edge colours
    dat_plt.MarkerEdgeColors = dat_plt.Colors ;
    % marker face colours
    c = cell(1, nDat) ;
    c(:) = {'none'} ;
    dat_plt.MarkerFaceColors = c ;
    % marker face alpha
    c = cell(1, nDat) ;
    c(:) = {1};
    dat_plt.MarkerAlpha = c ;

    % aspect ratio of the data
    dat_plt.ASR = 1/aspect_ratio ;
    % width of the figure
    dat_plt.figWidth_cm = baseWdith *plotMag ; % [cm]
    % fontsize
    dat_plt.AxesFontSize = 8 *plotMag ;
    % interpreter
    dat_plt.interpreter = 'tex' ;
    % xlabel
    dat_plt.XLabel = '{\itt}{\fontsize{8} }(s)' ;
    % ylabel
    dat_plt.YLabel = '{\it\delta}_{a}{\fontsize{8} }(nm)' ;
    % x limit range
    a = dat_plt.dat{1} ;
    b = dat_plt.dat{end} ;
    dat_plt.XLim = [0.01, 1.1*max(a(:,1))] ;
    % y limit range
    dat_plt.YLim = [0, 300] ; %[0, 1.1*max(b(:,2))] ;
    % number of x-axis labels [negative, positive]
    dat_plt.XLabNum = [1, 3] ;
    dat_plt.XTickValues = [1e-2, 0.1, 1, 10, max(a(:,1)) ] ; %1e-4, 1e-3, 
    dat_plt.XTickLabels = {'0.01', '0.1', '1', '10', num2str( max(a(:,1)), '%.0f' )} ;
    % number of y-axis labels [negative, positive]
    dat_plt.YLabNum = [1, 5] ;
    % x-aixs labels' format
    dat_plt.XLabForm = '%.0f' ;
    % y-aixs labels' format
    dat_plt.YLabForm = '%.0f' ;
    % x-axis scale type
    dat_plt.XScale = 'log' ;
    % y-axis scale type
    dat_plt.YScale = 'linear' ;
    % axis direction
    dat_plt.AxDir = 'left' ;
    % type of axis
    dat_plt.axType = 'single' ;

    dat_plt.splin = 'no' ;
    dat_plt.er_num = 3 ;

    % plot the data
    RA_Plot_1D(dat_plt) ;
    grid on

    % legend
    leg_tex = cell(nW, 1) ;
    for i = 1 : nW

        f = 1.0 + 0.05 * (i-1) ;
        leg_tex{i} = num2str(f, '%.2f') ;

    end
    leg_tex = leg_tex(end:-1:1) ;
    legend(leg_tex, 'fontsize', dat_plt.AxesFontSize, 'interpreter', 'tex', 'location', 'northwest', 'NumColumns', 1 ) ;

    % add text
    % text(0.5, 120, ['10^3\langle\beta\rangle{\fontsize{3} }={\fontsize{3} }', num2str(beta*1e3, '%.0f'), '\pm', num2str(sig_beta*1e3, '%.0f') ], 'interpreter', 'tex', 'fontsize', dat_plt.AxesFontSize ) ;
    % text(1.2*min(DRange), 1.1*min(URange), ['{\itD}_d{\fontsize{3} }={\fontsize{3} }', num2str(DRange(1), '%.0f'), '-', num2str(DRange(2), '%.0f'), '{\fontsize{3} }', char(181), 'm' ], 'interpreter', 'tex', 'fontsize', dat_plt.AxesFontSize ) ;

    %----------------------------------------------------------------------
    % Export the figure as a vectorised PDF
    %----------------------------------------------------------------------
    exportgraphics( gca, 'PDMS-20cst-SiO_1.00_SiW_varies.svg', 'ContentType', 'vector') ;
end
clear dat_plt
close

%% Figure C: Hexadecane | oil saturation varies

is_plt = 0 ;
if isequal( is_plt, 1 )
    %----------------------------------------------------------------------
    % plot properties
    %----------------------------------------------------------------------
    % number of data sets
    nDat = 7 ;
    % dataset collection
    dat_plt.dat = delta_nm_Hex_SoVaries(end:-1:1) ;
    % symbolsOS_PU
    dat_plt.Symbols = { '-', '--', '-.', '-', '--', '-.', '-' } ;
    % line colours
    dat_plt.Colors = {RA_getColorHex('Teal'), RA_getColorHex('DarkCyan'), RA_getColorHex('CadetBlue'), ...
        RA_getColorHex('SteelBlue'), RA_getColorHex('LightSeaGreen'),...
        RA_getColorHex('MediumAquamarine'), RA_getColorHex('DarkSlateGray')} ;
    % linewidths
    c = cell(1, nDat) ;
    c(:) = {0.75*plotMag} ;
    dat_plt.LineWidths = c ;
    % line alpha
    c = cell(1, nDat) ;
    c(:) = {1} ;
    dat_plt.LineAlpha = c ;

    % marker sizes
    c = cell(1, nDat) ;
    c(:) = {5*plotMag} ;
    dat_plt.MarkerSizes = c ;
    % edge colours
    dat_plt.MarkerEdgeColors = dat_plt.Colors ;
    % marker face colours
    c = cell(1, nDat) ;
    c(:) = {'none'} ;
    dat_plt.MarkerFaceColors = c ;
    % marker face alpha
    c = cell(1, nDat) ;
    c(:) = {1};
    dat_plt.MarkerAlpha = c ;

    % aspect ratio of the data
    dat_plt.ASR = 1/aspect_ratio ;
    % width of the figure
    dat_plt.figWidth_cm = baseWdith *plotMag ; % [cm]
    % fontsize
    dat_plt.AxesFontSize = 8 *plotMag ;
    % interpreter
    dat_plt.interpreter = 'tex' ;
    % xlabel
    dat_plt.XLabel = '{\itt}{\fontsize{8} }(s)' ;
    % ylabel
    dat_plt.YLabel = '{\it\delta}_{a}{\fontsize{8} }(nm)' ;
    % x limit range
    a = dat_plt.dat{1} ;
    b = dat_plt.dat{end} ;
    dat_plt.XLim = [1e-4, 1.1*max(a(:,1))] ;
    % y limit range
    dat_plt.YLim = [0, 200] ; %[0.9*min(b(:,2)), 1.1*max(b(:,2))] ;
    % number of x-axis labels [negative, positive]
    dat_plt.XLabNum = [1, 3] ;
    dat_plt.XTickValues = [1e-4, 1e-3, 1e-2, 0.1, 1, 10] ;
    dat_plt.XTickLabels = {'10^{-4}', '10^{-3}', '10^{-2}', '0.1', '1', '10'} ;
    % number of y-axis labels [negative, positive]
    dat_plt.YLabNum = [1, 5] ;
    % x-aixs labels' format
    dat_plt.XLabForm = '%.0f' ;
    % y-aixs labels' format
    dat_plt.YLabForm = '%.0f' ;
    % x-axis scale type
    dat_plt.XScale = 'log' ;
    % y-axis scale type
    dat_plt.YScale = 'linear' ;
    % axis direction
    dat_plt.AxDir = 'left' ;
    % type of axis
    dat_plt.axType = 'single' ;

    dat_plt.splin = 'no' ;
    dat_plt.er_num = 3 ;

    % plot the data
    RA_Plot_1D(dat_plt) ;
    grid on

    % legend
    leg_tex = cell(nW, 1) ;
    for i = 1 : nW

        f = 1.0 + 0.05 * (i-1) ;
        leg_tex{i} = num2str(f, '%.2f') ;

    end
    leg_tex = leg_tex(end:-1:1) ;
    legend(leg_tex, 'fontsize', dat_plt.AxesFontSize, 'interpreter', 'tex', 'location', 'northwest', 'NumColumns', 1 ) ;

    % add text
    % text(0.5, 120, ['10^3\langle\beta\rangle{\fontsize{3} }={\fontsize{3} }', num2str(beta*1e3, '%.0f'), '\pm', num2str(sig_beta*1e3, '%.0f') ], 'interpreter', 'tex', 'fontsize', dat_plt.AxesFontSize ) ;
    % text(1.2*min(DRange), 1.1*min(URange), ['{\itD}_d{\fontsize{3} }={\fontsize{3} }', num2str(DRange(1), '%.0f'), '-', num2str(DRange(2), '%.0f'), '{\fontsize{3} }', char(181), 'm' ], 'interpreter', 'tex', 'fontsize', dat_plt.AxesFontSize ) ;

    %----------------------------------------------------------------------
    % Export the figure as a vectorised PDF
    %----------------------------------------------------------------------
    exportgraphics( gca, 'Hexadecane-SiO_varies_SiW_1.00.svg', 'ContentType', 'vector') ;
end
clear dat_plt
close

%% Figure D: 20 cst PDMS | oil saturation varies

is_plt = 0 ;
if isequal( is_plt, 1 )
    %----------------------------------------------------------------------
    % plot properties
    %----------------------------------------------------------------------
    % number of data sets
    nDat = 7 ;
    % dataset collection
    dat_plt.dat = delta_nm_20cst_SoVaries(end:-1:1) ;
    % symbolsOS_PU
    dat_plt.Symbols = { '-', '--', '-.', '-', '--', '-.', '-' } ;
    % line colours
    dat_plt.Colors = {RA_getColorHex('DarkRed'), RA_getColorHex('FireBrick'), RA_getColorHex('Crimson'), ...
        RA_getColorHex('IndianRed'), RA_getColorHex('Tomato'),...
        RA_getColorHex('OrangeRed'), RA_getColorHex('DarkOrange')} ;
    % linewidths
    c = cell(1, nDat) ;
    c(:) = {0.75*plotMag} ;
    dat_plt.LineWidths = c ;
    % line alpha
    c = cell(1, nDat) ;
    c(:) = {1} ;
    dat_plt.LineAlpha = c ;

    % marker sizes
    c = cell(1, nDat) ;
    c(:) = {5*plotMag} ;
    dat_plt.MarkerSizes = c ;
    % edge colours
    dat_plt.MarkerEdgeColors = dat_plt.Colors ;
    % marker face colours
    c = cell(1, nDat) ;
    c(:) = {'none'} ;
    dat_plt.MarkerFaceColors = c ;
    % marker face alpha
    c = cell(1, nDat) ;
    c(:) = {1};
    dat_plt.MarkerAlpha = c ;

    % aspect ratio of the data
    dat_plt.ASR = 1/aspect_ratio ;
    % width of the figure
    dat_plt.figWidth_cm = baseWdith *plotMag ; % [cm]
    % fontsize
    dat_plt.AxesFontSize = 8 *plotMag ;
    % interpreter
    dat_plt.interpreter = 'tex' ;
    % xlabel
    dat_plt.XLabel = '{\itt}{\fontsize{8} }(s)' ;
    % ylabel
    dat_plt.YLabel = '{\it\delta}_{a}{\fontsize{8} }(nm)' ;
    % x limit range
    a = dat_plt.dat{1} ;
    b = dat_plt.dat{end} ;
    dat_plt.XLim = [1e-4, 1.1*max(a(:,1))] ;
    % y limit range
    dat_plt.YLim = [0, 400] ; %[0.9*min(b(:,2)), 1.1*max(b(:,2))] ;
    % number of x-axis labels [negative, positive]
    dat_plt.XLabNum = [1, 3] ;
    dat_plt.XTickValues = [1e-4, 1e-3, 1e-2, 0.1, 1, 10, max(a(:,1)) ] ;
    dat_plt.XTickLabels = {'10^{-4}', '10^{-3}', '10^{-2}', '0.1', '1', '10', num2str( max(a(:,1)), '%.0f' )} ;
    % number of y-axis labels [negative, positive]
    dat_plt.YLabNum = [1, 5] ;
    % x-aixs labels' format
    dat_plt.XLabForm = '%.0f' ;
    % y-aixs labels' format
    dat_plt.YLabForm = '%.0f' ;
    % x-axis scale type
    dat_plt.XScale = 'log' ;
    % y-axis scale type
    dat_plt.YScale = 'linear' ;
    % axis direction
    dat_plt.AxDir = 'left' ;
    % type of axis
    dat_plt.axType = 'single' ;

    dat_plt.splin = 'no' ;
    dat_plt.er_num = 3 ;

    % plot the data
    RA_Plot_1D(dat_plt) ;
    grid on

    % legend
    leg_tex = cell(nW, 1) ;
    for i = 1 : nW

        f = 1.0 + 0.05 * (i-1) ;
        leg_tex{i} = num2str(f, '%.2f') ;

    end
    leg_tex = leg_tex(end:-1:1) ;
    legend(leg_tex, 'fontsize', dat_plt.AxesFontSize, 'interpreter', 'tex', 'location', 'northwest', 'NumColumns', 1 ) ;

    % add text
    % text(0.5, 120, ['10^3\langle\beta\rangle{\fontsize{3} }={\fontsize{3} }', num2str(beta*1e3, '%.0f'), '\pm', num2str(sig_beta*1e3, '%.0f') ], 'interpreter', 'tex', 'fontsize', dat_plt.AxesFontSize ) ;
    % text(1.2*min(DRange), 1.1*min(URange), ['{\itD}_d{\fontsize{3} }={\fontsize{3} }', num2str(DRange(1), '%.0f'), '-', num2str(DRange(2), '%.0f'), '{\fontsize{3} }', char(181), 'm' ], 'interpreter', 'tex', 'fontsize', dat_plt.AxesFontSize ) ;

    %----------------------------------------------------------------------
    % Export the figure as a vectorised PDF
    %----------------------------------------------------------------------
    exportgraphics( gca, 'PDMS-20cst-SiO_varies_SiW_1.00.svg', 'ContentType', 'vector') ;
end
close
clear dat_plt

%% Figure E: Hexadecane | water saturation & droplet diameter vary

is_plt = 0 ;
if isequal( is_plt, 1 )
    %----------------------------------------------------------------------
    % arrange the data
    %----------------------------------------------------------------------
    Dd = (150:10:400)' ; % [um]
    dat = cell(nW, 1) ;
    for i = 1 : nW
        dat{i,1} = [Dd, delta_Vis_Hex_SwVaries(:,i)] ;
    end
    dat = dat(end:-1:1) ;
    %----------------------------------------------------------------------
    % plot properties
    %----------------------------------------------------------------------
    % number of data sets
    nDat = 7 ;
    % dataset collection
    dat_plt.dat = dat ;
    % symbolsOS_PU
    dat_plt.Symbols = { '-', '--', '-.', ':', '--', '-.', '-' } ;
    % line colours
    dat_plt.Colors = {RA_getColorHex('Teal'), RA_getColorHex('DarkCyan'), RA_getColorHex('CadetBlue'), ...
        RA_getColorHex('SteelBlue'), RA_getColorHex('LightSeaGreen'),...
        RA_getColorHex('MediumAquamarine'), RA_getColorHex('DarkSlateGray')} ;
    % linewidths
    c = cell(1, nDat) ;
    c(:) = {0.75*plotMag} ;
    dat_plt.LineWidths = c ;
    % line alpha
    c = cell(1, nDat) ;
    c(:) = {1} ;
    dat_plt.LineAlpha = c ;

    % marker sizes
    c = cell(1, nDat) ;
    c(:) = {5*plotMag} ;
    dat_plt.MarkerSizes = c ;
    % edge colours
    dat_plt.MarkerEdgeColors = dat_plt.Colors ;
    % marker face colours
    c = cell(1, nDat) ;
    c(:) = {'none'} ;
    dat_plt.MarkerFaceColors = c ;
    % marker face alpha
    c = cell(1, nDat) ;
    c(:) = {1};
    dat_plt.MarkerAlpha = c ;

    % aspect ratio of the data
    dat_plt.ASR = 1/aspect_ratio ;
    % width of the figure
    dat_plt.figWidth_cm = baseWdith *plotMag ; % [cm]
    % fontsize
    dat_plt.AxesFontSize = 8 *plotMag ;
    % interpreter
    dat_plt.interpreter = 'tex' ;
    % xlabel
    dat_plt.XLabel = ['{\itD}_{d}{\fontsize{8} }(', char(181), 'm)'] ;
    % ylabel
    dat_plt.YLabel = '{\it\delta}_{a,\infty}{\fontsize{8} }(nm)' ;
    % x limit range
    a = dat_plt.dat{1} ;
    b = dat_plt.dat{end} ;
    dat_plt.XLim = [150, 400] ;
    % y limit range
    dat_plt.YLim = [0, 200] ; %[0.9*min(a(:,2)), 1.1*max(b(:,2))] ;
    % number of x-axis labels [negative, positive]
    dat_plt.XLabNum = [1, 4] ;
    % dat_plt.XTickValues = [1e-4, 1e-3, 1e-2, 0.1, 1, 10, max(a(:,2)) ] ;
    % dat_plt.XTickLabels = {'10^{-4}', '10^{-3}', '10^{-2}', '0.1', '10', num2str( max(a(:,1)), '%.0f' )} ;
    % number of y-axis labels [negative, positive]
    dat_plt.YLabNum = [1, 4] ;
    % x-aixs labels' format
    dat_plt.XLabForm = '%.0f' ;
    % y-aixs labels' format
    dat_plt.YLabForm = '%.0f' ;
    % x-axis scale type
    dat_plt.XScale = 'linear' ;
    % y-axis scale type
    dat_plt.YScale = 'linear' ;
    % axis direction
    dat_plt.AxDir = 'left' ;
    % type of axis
    dat_plt.axType = 'single' ;

    dat_plt.splin = 'no' ;
    dat_plt.er_num = 3 ;

    % plot the data
    RA_Plot_1D(dat_plt) ;
    grid on
    % grid minor

    % legend
    leg_tex = cell(nW, 1) ;
    for i = 1 : nW

        f = 1 + 0.05 * (i-1) ;
        leg_tex{i} = num2str(f, '%.2f') ;

    end
    leg_tex = leg_tex(end:-1:1) ;
    legend(leg_tex, 'fontsize', dat_plt.AxesFontSize, 'interpreter', 'tex', 'location', 'northeast', 'NumColumns', 2, 'Box','off' ) ;

    % add text
    % text(0.5, 120, ['10^3\langle\beta\rangle{\fontsize{3} }={\fontsize{3} }', num2str(beta*1e3, '%.0f'), '\pm', num2str(sig_beta*1e3, '%.0f') ], 'interpreter', 'tex', 'fontsize', dat_plt.AxesFontSize ) ;
    % text(1.2*min(DRange), 1.1*min(URange), ['{\itD}_d{\fontsize{3} }={\fontsize{3} }', num2str(DRange(1), '%.0f'), '-', num2str(DRange(2), '%.0f'), '{\fontsize{3} }', char(181), 'm' ], 'interpreter', 'tex', 'fontsize', dat_plt.AxesFontSize ) ;

    %----------------------------------------------------------------------
    % Export the figure as a vectorised PDF
    %----------------------------------------------------------------------
    exportgraphics( gca, 'Hexadecane-SiO_1.00_SiW_&_Diameter_varies.svg', 'ContentType', 'vector') ;
end
clear dat_plt
close

%% Figure F: 20 cst PDMS | water saturation & droplet diameter vary

is_plt = 0 ;
if isequal( is_plt, 1 )
    %----------------------------------------------------------------------
    % arrange the data
    %----------------------------------------------------------------------
    Dd = (150:10:400)' ; % [um]
    dat = cell(nW, 1) ;
    for i = 1 : nW
        dat{i,1} = [Dd, delta_Vis_20cst_SwVaries(:,i)] ;
    end
    dat = dat(end:-1:1) ;
    %----------------------------------------------------------------------
    % plot properties
    %----------------------------------------------------------------------
    % number of data sets
    nDat = 7 ;
    % dataset collection
    dat_plt.dat = dat ;
    % symbolsOS_PU
    dat_plt.Symbols = { '-', '--', '-.', ':', '--', '-.', '-' } ;
    % line colours
    dat_plt.Colors = {RA_getColorHex('FireBrick'), RA_getColorHex('IndianRed'), RA_getColorHex('Tomato'), ...
        RA_getColorHex('Coral'), RA_getColorHex('Chocolate'),...
        RA_getColorHex('Sienna'), RA_getColorHex('DarkOrange')} ;
    % linewidths
    c = cell(1, nDat) ;
    c(:) = {0.75*plotMag} ;
    dat_plt.LineWidths = c ;
    % line alpha
    c = cell(1, nDat) ;
    c(:) = {1} ;
    dat_plt.LineAlpha = c ;

    % marker sizes
    c = cell(1, nDat) ;
    c(:) = {5*plotMag} ;
    dat_plt.MarkerSizes = c ;
    % edge colours
    dat_plt.MarkerEdgeColors = dat_plt.Colors ;
    % marker face colours
    c = cell(1, nDat) ;
    c(:) = {'none'} ;
    dat_plt.MarkerFaceColors = c ;
    % marker face alpha
    c = cell(1, nDat) ;
    c(:) = {1};
    dat_plt.MarkerAlpha = c ;

    % aspect ratio of the data
    dat_plt.ASR = 1/aspect_ratio ;
    % width of the figure
    dat_plt.figWidth_cm = baseWdith *plotMag ; % [cm]
    % fontsize
    dat_plt.AxesFontSize = 8 *plotMag ;
    % interpreter
    dat_plt.interpreter = 'tex' ;
    % xlabel
    dat_plt.XLabel = ['{\itD}_{d}{\fontsize{8} }(', char(181), 'm)'] ;
    % ylabel
    dat_plt.YLabel = '{\it\delta}_{a,\infty}{\fontsize{8} }(nm)' ;
    % x limit range
    a = dat_plt.dat{1} ;
    b = dat_plt.dat{end} ;
    dat_plt.XLim = [150, 400] ;
    % y limit range
    dat_plt.YLim = [0, 400] ; %[0.9*min(a(:,2)), 1.1*max(b(:,2))] ;
    % number of x-axis labels [negative, positive]
    dat_plt.XLabNum = [1, 4] ;
    % dat_plt.XTickValues = [1e-4, 1e-3, 1e-2, 0.1, 1, 10, max(a(:,2)) ] ;
    % dat_plt.XTickLabels = {'10^{-4}', '10^{-3}', '10^{-2}', '0.1', '10', num2str( max(a(:,1)), '%.0f' )} ;
    % number of y-axis labels [negative, positive]
    dat_plt.YLabNum = [1, 4] ;
    % x-aixs labels' format
    dat_plt.XLabForm = '%.0f' ;
    % y-aixs labels' format
    dat_plt.YLabForm = '%.0f' ;
    % x-axis scale type
    dat_plt.XScale = 'linear' ;
    % y-axis scale type
    dat_plt.YScale = 'linear' ;
    % axis direction
    dat_plt.AxDir = 'left' ;
    % type of axis
    dat_plt.axType = 'single' ;

    dat_plt.splin = 'no' ;
    dat_plt.er_num = 3 ;

    % plot the data
    RA_Plot_1D(dat_plt) ;
    grid on
    % grid minor

    % legend
    leg_tex = cell(nW, 1) ;
    for i = 1 : nW

        f = 1.00 + 0.05 * (i-1) ;
        leg_tex{i} = num2str(f, '%.2f') ;

    end
    leg_tex = leg_tex(end:-1:1) ;
    legend(leg_tex, 'fontsize', dat_plt.AxesFontSize, 'interpreter', 'tex', 'location', 'northeast', 'NumColumns', 2, 'box', 'off' ) ;

    % add text
    % text(0.5, 120, ['10^3\langle\beta\rangle{\fontsize{3} }={\fontsize{3} }', num2str(beta*1e3, '%.0f'), '\pm', num2str(sig_beta*1e3, '%.0f') ], 'interpreter', 'tex', 'fontsize', dat_plt.AxesFontSize ) ;
    % text(1.2*min(DRange), 1.1*min(URange), ['{\itD}_d{\fontsize{3} }={\fontsize{3} }', num2str(DRange(1), '%.0f'), '-', num2str(DRange(2), '%.0f'), '{\fontsize{3} }', char(181), 'm' ], 'interpreter', 'tex', 'fontsize', dat_plt.AxesFontSize ) ;

    %----------------------------------------------------------------------
    % Export the figure as a vectorised PDF
    %----------------------------------------------------------------------
    exportgraphics( gca, 'PDMS-20cst-SiO_1.00_SiW_&_Diameter_varies.svg', 'ContentType', 'vector') ;
end
clear dat_plt


