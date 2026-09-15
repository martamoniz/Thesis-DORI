%% ======================================================================
%  POSTER VERSION of the Co-57 Compton edge, on the channel axis
%
%  Bottom left picture of the three-figure panel in RESULTS. Same size,
%  same fonts and same line weights as am241_sub.m, so the pair reads as
%  one figure.
%
%  Same construction as the thesis figure: the zoom window around the edge
%  is drawn as connected markers, and the centroid is marked with a yellow
%  square carrying its horizontal uncertainty. Nothing is fitted here; the
%  centroids are the ones read off the five acquisitions.
%
%  Built at the size it occupies on the poster, so LaTeX includes it at
%  width=\subwid and applies no scaling.
% ======================================================================

W_CM = 11.6;   H_CM = 10.8;       % final size on the poster, centimetres

FS_TICK = 20;                     % tick labels
FS_LAB  = 23;                     % axis labels
FS_LEG  = 17;                     % legend
FS_TTL  = 23;                     % title

LW_AXES = 2.0;
LW_DATA = 2.0;                    % spectrum
LW_ERR  = 3.5;                    % centroid error bar
MS_PT   = 5;                      % data markers
MS_MU   = 150;                    % centroid marker, scatter units

col_yellow = [204 153 0]/255;     % the yellow the Co-57 point uses in the
dark       = [0 0 0];             % calibration figure

DOTS_LO = 34;   ZOOM_HI = 45;     % zoom window around the edge, channels

%% ------------------------------------------------------------- data
srcPath = 'C:\Users\marta\OneDrive - Universidade de Aveiro\DORI_25-26\RESULTS\spectrum\1_junho\Sources\';

bg       = readtable([srcPath 'background.csv']);
bgCounts = bg.counts;

data0     = readtable([srcPath 'Co-57.csv']);
bins      = data0.bin;
corrected = data0.counts - bgCounts;
corrected(corrected < 0) = 0;

% Centroids read off the five acquisitions, as in the calibration script.
centroids = [39.0; 37.0; 38.0; 39.0; 39.0];
mu_bar    = median(centroids);
muErr     = max(abs(centroids - mu_bar));

maskD = bins >= DOTS_LO & bins <= ZOOM_HI;
xData = bins(maskD);
yData = corrected(maskD);

xLeft  = DOTS_LO - 0.5;
xRight = ZOOM_HI + 0.5;
yhi_ax = max(yData)*1.1;
ylo_ax = min(yData(yData > 0))*0.9;

step    = max(round((yhi_ax - ylo_ax)/3/100)*100, 50);
y_major = (ceil(ylo_ax/step)*step) : step : (floor(yhi_ax/step)*step);
x_major = (ceil(xLeft/5)*5) : 5 : (floor(xRight/5)*5);
x_minor = ceil(xLeft) : 1 : floor(xRight);

%% ------------------------------------------------------------- figure
figP = figure('Color','w', 'Units','centimeters', ...
              'Position',[2 2 W_CM H_CM], ...
              'PaperUnits','centimeters', 'PaperPosition',[0 0 W_CM H_CM]);
axP = axes(figP);

h_dat = plot(axP, xData, yData, 'o-', 'MarkerSize', MS_PT, ...
             'MarkerEdgeColor', dark, 'MarkerFaceColor', dark, ...
             'Color', dark, 'LineWidth', LW_DATA, ...
             'DisplayName', '$^{57}$Co spectrum');
hold(axP,'on');

yBar   = max(yData);
mu_lbl = sprintf('$C_E = %.0f \\pm %.0f$ (a.u.)', mu_bar, muErr);

errorbar(axP, mu_bar, yBar, muErr, 'horizontal', ...
         'Color', col_yellow, 'LineWidth', LW_ERR, 'CapSize', 10, ...
         'LineStyle','none', 'Marker','none', 'HandleVisibility','off');

h_mu = scatter(axP, mu_bar, yBar, MS_MU, 's', ...
               'MarkerFaceColor', col_yellow, 'MarkerEdgeColor', col_yellow, ...
               'MarkerFaceAlpha', 0.6, 'LineWidth', 1.5, ...
               'DisplayName', mu_lbl);

hold(axP,'off');

set(axP, ...
    'TickLabelInterpreter','latex', ...
    'FontSize',   FS_TICK, ...
    'LineWidth',  LW_AXES, ...
    'TickDir',    'in', ...
    'TickLength', [0.018 0.010], ...
    'XLim', [xLeft xRight], 'XTick', x_major, 'XTickMode','manual', ...
    'YLim', [ylo_ax yhi_ax], 'YTick', y_major, 'YTickMode','manual', ...
    'XMinorTick','on', 'YMinorTick','off', ...
    'Box','on');
axP.XAxis.MinorTickValues = x_minor;
grid(axP,'off');                  % no grid on the poster version

xlabel(axP, 'Channel (a.u.)', 'Interpreter','latex', 'FontSize', FS_LAB);
ylabel(axP, 'Counts (a.u.)',  'Interpreter','latex', 'FontSize', FS_LAB);
title(axP,  '\textbf{$^{57}$Co}', 'Interpreter','latex', ...
      'FontSize', FS_TTL, 'Color', dark);

% The axes are placed before the legend is created, so the legend anchors
% itself inside the final plot box instead of the default one.
set(axP, 'Units','normalized', 'Position', [0.215 0.195 0.755 0.680]);

legP = legend(axP, [h_dat h_mu], 'Interpreter','latex', ...
              'Location','southwest', 'Box','on', 'FontSize', FS_LEG);
legP.LineWidth     = 1.2;
legP.ItemTokenSize = [26 18];

% Widen the legend after the LaTeX text has been rendered, keeping the left
% edge where it was.
drawnow;
legP.Position(3) = legP.Position(3)*1.20;

fprintf('Co-57: centroid = %.0f +/- %.0f channels\n', mu_bar, muErr);

exportgraphics(figP, 'co57_sub.pdf', ...
    'ContentType','vector', 'BackgroundColor','white');
fprintf('Exported: co57_sub.pdf  (%.1f x %.1f cm)\n', W_CM, H_CM);
