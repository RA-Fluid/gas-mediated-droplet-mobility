function RA_Plot_1D(Dat)
% This function plots the 2D results, whether it is an image, scalar, or a
% velocity field.
% What to be displayed is fed into the function by the user
%
% Developed by Reza Azadi, UoAlberta, (2021)
%
% Plt_Type: Im (Image), Scal (Scalar field), Vel (Velocity, 2D)
%
% Dat.X                 : x-positions in [mm]
% Dat.Y                 : y-positions in [mm]
% Dat.ASR               : Data aspect ratio
% Dat.AxesFontSize      : Axes font size, e.g., 40
% Dat.XLabel            : X-label string
% Dat.YLabel            : Y-label string

%% Initialise

% open a figure
if isequal( Dat.AxDir, 'left' )
    f = figure(1);
    set(f, 'Units', 'centimeters');
    set(f, 'Position', [1, 2, Dat.figWidth_cm, Dat.figWidth_cm * Dat.ASR] );
    % fig = figure('unit', 'normalised', 'outerposition', [0,0,0.75,1], 'visible', 'on') ;

    left_color = [0 0 0];
    right_color = [0 0 0];
    set(f,'defaultAxesColorOrder', [left_color; right_color] );

    % set background color to white
    set(gcf, 'color', 'w' ) ;
    hold on ;
end

if isequal( Dat.axType, 'dual')
    dir = Dat.AxDir ;
    yyaxis(dir) ;
end

%% get the axes information

ax = gca ;

% Set the color of the axes lines and tick labels to black
ax.XColor = 'k';
ax.YColor = 'k';

% Set the color of the text (title and labels) to black
ax.Title.Color = 'k';
ax.XLabel.Color = 'k';
ax.YLabel.Color = 'k';

%the data limits
XLim = Dat.XLim ;
YLim = Dat.YLim ;

% read the data: each array of the cell is one set of data, with the first
% column the x and second column the y data.
dat = Dat.dat ;
sym = Dat.Symbols ; % symbols
col = Dat.Colors ; % line colors
linwid = Dat.LineWidths ; % line widths
marksize = Dat.MarkerSizes ; % marker sizes
markedgcol = Dat.MarkerEdgeColors ; % marker edge colors
markfacol = Dat.MarkerFaceColors ; % marker face colors

%% plot

for i = 1 : length(dat)

    % extract the x & y data
    x = dat{i}(:,1) ;
    y = dat{i}(:,2) ;

    % data with no error bars
    if size( dat{i}, 2 ) == 2

        dx = [] ;
        dy = [] ;
        % xs = linspace(x(1), x(end), 10*length(x)) ;
        % ys = pchip(x,y, xs) ;
        % x = xs ;
        % y = ys ;
    % data with only y-errorbars
    elseif size( dat{i}, 2 ) == 3
        dx = nan(size(dy)) ;
        dy = dat{i}(:,3) ;
    % data with x-& y-errorbars
    else
        dx = dat{i}(:, 3) ;
        dy = dat{i}(:, 4) ;        
    end

    if isequal( markfacol{i}, 'none')
        sym_col_alpha = markfacol{i} ;
    else
        sym_col_alpha = RA_getAlphaColor( markfacol{i}, Dat.MarkerAlpha{i}, [1 1 1] ) ;
    end
    lin_col_alpha = RA_getAlphaColor( col{i}, Dat.LineAlpha{i}, [1 1 1] ) ;

    % errorbar( x, y, dy, dy, dx, dx, sym{i}, 'Color', lin_col_alpha, 'linewidth', linwid{i}, ...
    %     'Markersize', marksize{i}, 'markeredgecolor', sym_col_alpha,...
    %     'markerfacecolor', markfacol{i} ) ;

    %----------------------------------------------------------------------
    % plot the errorbars first and the symbols over it
    %----------------------------------------------------------------------
    hold on;

    % Draw the Symbols (Foreground)
    % We use standard 'plot' to put the markers on top of the bars.
    plot( x, y, sym{i}, ...
        'Color', lin_col_alpha, ...          % Keep connecting line color same (if sym has line)
        'LineWidth', linwid{i}, ...
        'MarkerSize', marksize{i}, ...
        'MarkerEdgeColor', lin_col_alpha, ...%markedgcol{i}
        'MarkerFaceColor', sym_col_alpha );
    
    % STEP 1: Draw the Error Bars (Background)
    % We use your existing errorbar command, but force 'Marker' to 'none'.
    % This draws ONLY the vertical/horizontal bars and the connecting line.
    if ~isequal(sym{i}, '-') && ~isequal(sym{i}, '--') && ~isequal(sym{i}, ':') && ~isequal(sym{i}, '-.')
        h_err = errorbar( x, y, dy, dy, dx, dx, sym{i}, ...
            'Color', lin_col_alpha, ...%markedgcol{i}
            'LineWidth', 0.5*linwid{i}, ...
            'Marker', 'none' ); %Turn off markers here
        
        % Prevent this layer from appearing in the legend
        set(get(get(h_err,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end

    % hold off;

end
% plot a spline
if isequal( Dat.splin, 'yes' )
    for i = 1 : length(dat)
        
        x = dat{i}(:,1) ;
        y = dat{i}(:,2) ;
        if i <= Dat.er_num

            xs = linspace(x(1), x(end), 5*length(x)) ;
            ys = spline(x,y, xs) ;
            x = xs ;
            y = ys ;

            sym_col_alpha = RA_getAlphaColor(col{i}, Dat.LineAlpha, [1 1 1]) ;

            plot( x, y, '--', 'Color', sym_col_alpha, 'linewidth', linwid{i}, ...
                'Markersize', marksize{i}, 'markeredgecolor', markedgcol{i},...
                'markerfacecolor', markfacol{i} ) ;
        end

    end
end

%% Set the global properties

% set axes properties
pbaspect([ 1 Dat.ASR 1] ) ;
% ax.FontName = 'Times New Roman' ;
ax.FontSize = Dat.AxesFontSize ;
ax.TickDir = 'out' ;
ax.TickLabelInterpreter = Dat.interpreter ;
ax.Box = 'on' ;
ax.LineWidth = 0.8*linwid{1}*0.25 ;
ax.XScale = Dat.XScale ;
ax.YScale = Dat.YScale ;

ax.XLabel.Interpreter = Dat.interpreter ;
ax.YLabel.Interpreter = Dat.interpreter ;
ax.XLabel.String = Dat.XLabel ;
ax.YLabel.String = Dat.YLabel ;
ax.XLim = XLim ;
ax.YLim = YLim ;

% set x-axis ticklabels
if isequal( Dat.XScale, 'linear')

    if XLim(1) < 0 && XLim(2) > 0
        nx = [ ( XLim(1) : -XLim(1)/Dat.XLabNum(1) : 0 ), ( XLim(2)/Dat.XLabNum(2) : XLim(2)/Dat.XLabNum(2) : XLim(2) ) ] ;
    elseif XLim(1) < 0 && XLim(2) < 0
        nx = XLim(1) : diff(XLim) / Dat.XLabNum(1) : XLim(2) ;
    elseif XLim(1) >= 0 && XLim(2) >= 0
        nx = XLim(1) : diff(XLim) / Dat.XLabNum(2) : XLim(2) ;
    end

    xlab = cell(size(nx)) ;
    for i = 1 : length(nx)
        if nx(i) == 0
            xlab(i) = arrayfun( @(x) sprintf( '%.0f', x), nx(i), 'un', 0) ;
        else
            xlab(i) = arrayfun( @(x) sprintf(  Dat.XLabForm, x), nx(i), 'un', 0) ;
        end

    end

    ax.XTick = nx ; % x-axis tick values
    ax.XTickLabel = xlab ; % x-axis ticklabels

elseif isequal( Dat.XScale, 'log')

    ax.XTick = Dat.XTickValues ; % x-axis tick values
    ax.XTickLabel = Dat.XTickLabels ; % x-axis ticklabels
    % ax.XAxis.Exponent = 3;
end

% set y-axis ticklabels
if isequal( Dat.YScale, 'linear')

    if YLim(1)< 0
        ny = [ ( YLim(1) : -YLim(1)/Dat.YLabNum(1) : 0 ), ( YLim(2)/Dat.YLabNum(2) : YLim(2)/Dat.YLabNum(2) : YLim(2) ) ] ;
    else
        ny = YLim(1): diff(YLim)/Dat.YLabNum(2) : YLim(2) ;
    end

    ylab = cell(size(ny)) ;
    for i = 1 : length(ny)
        if ny(i) == 0
            ylab(i) = arrayfun( @(x) sprintf( '%.0f', x), ny(i), 'un', 0) ;
        else
            if isequal( Dat.YScale, 'linear')
                ylab(i) = arrayfun( @(x) sprintf( Dat.YLabForm, x), ny(i), 'un', 0) ;
            elseif isequal( Dat.YScale, 'log')
                ylab(i) = arrayfun( @(x) sprintf( Dat.YLabForm, x), ny(i), 'un', 0) ;
            end
        end

    end
    ax.YTick = ny ; % x-axis tick values
    ax.YTickLabel = ylab ; % x-axis ticklabels

elseif isequal( Dat.YScale, 'log')
    ax.YTick = Dat.YTickValues ; % x-axis tick values
    ax.YTickLabel = Dat.YTickLabels ; % x-axis ticklabels
    % ax.YAxis.Exponent = 3 ;
end