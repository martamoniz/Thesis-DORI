%% ======================================================================
%  POSTER VERSION of the Co-57 Compton edge spectrum, on the channel axis
%
%  Bottom left picture of the three-figure panel in RESULTS. Same size,
%  same fonts and same line weights as am241_sub.m, so the pair reads as
%  one figure.
%
%  TWO LINES TO CHECK, both marked EDIT below: the name of the .mat file
%  that holds the Co-57 fit, and the model the edge was fitted with. The
%  rest of the script is styling and does not depend on either. If the Co-57
%  fit is stored the same way as the Am-241 one (a fitResults struct with
%  bins, corrected, mu, fwhm and pHat), nothing needs to change.
% ======================================================================

W_CM = 11.6;   H_CM = 10.0;       % final size on the poster, centimetres

FS_TICK = 20;                     % tick labels
FS_LAB  = 23;                     % axis labels
FS_LEG  = 17;                     % legend
FS_TTL  = 23;                     % title

LW_AXES = 2.0;
LW_FIT  = 4;                      % fit curve
MS_PT   = 5;                      % data markers

col_yellow = [204 153 0]/255;     % same yellow the Co-57 point uses in the
dark       = [0 0 0];             % calibration figure

i = 1;                            % which acquisition to plot

%% ------------------------------------------------------------- data
srcPath = 'C:\Users\marta\OneDrive - Universidade de Aveiro\DORI_25-26\RESULTS\spectrum\1_junho\Sources\';

load([srcPath 'fitResults_Co57.mat'], 'fitResults');            % EDIT: file name

% EDIT: the model the Compton edge was fitted with. This is the Gaussian
% used for Am-241; swap it for the edge model if the Co-57 fit used one.
edge_model = @(pG, x) pG(1).*exp(-(x - pG(2)).^2 ./ (2*pG(3).^2)) + pG(4) + pG(5).*x;

bins      = double(fitResults(i).bins);
corrected = double(fitResults(i).corrected);
corrected(corrected < 1) = 1;                 % same floor as the thesis figure
pHat      = fitResults(i).pHat;

X_LO = min(bins);   X_HI = max(bins);         % channel window drawn
xCurve = linspace(double(fitResults(i).xFitLo), ...
                  double(fitResults(i).xFitHi), 500);
yCurve = edge_model(pHat, xCurve);
ytop   = 3*max(yCurve);

%% ------------------------------------------------------------- figure
figP = figure('Color','w', 'Units','centimeters', ...
              'Position',[2 2 W_CM H_CM], ...
              'PaperUnits','centimeters', 'PaperPosition',[0 0 W_CM H_CM]);
axP = axes(figP);
hold(axP,'on');

plot(axP, bins, corrected, 'o', 'MarkerSize', MS_PT, ...
     'MarkerEdgeColor', dark, 'MarkerFaceColor', dark, ...
     'HandleVisibility','off');

h_fit = plot(axP, xCurve, yCurve, '-', 'Color', col_yellow, 'LineWidth', LW_FIT);

hold(axP,'off');

set(axP, ...
    'TickLabelInterpreter','latex', ...
    'FontSize',   FS_TICK, ...
    'LineWidth',  LW_AXES, ...
    'TickDir',    'in', ...
    'TickLength', [0.018 0.010], ...
    'XLim', [X_LO X_HI], ...
    'YLim', [0 ytop], ...
    'XMinorTick','on', 'YMinorTick','off', ...
    'Box','on');
grid(axP,'off');                  % no grid on the poster version

xlabel(axP, 'Channel (a.u.)', 'Interpreter','latex', 'FontSize', FS_LAB);
ylabel(axP, 'Counts (a.u.)',  'Interpreter','latex', 'FontSize', FS_LAB);
title(axP,  '\textbf{$^{57}$Co}', 'Interpreter','latex', ...
      'FontSize', FS_TTL, 'Color', dark);

set(axP, 'Units','normalized', 'Position', [0.215 0.195 0.755 0.680]);

legP = legend(axP, h_fit, 'Compton Edge Fit', ...
              'Interpreter','latex', 'Location','northeast', ...
              'Box','on', 'FontSize', FS_LEG);
legP.LineWidth     = 1.2;
legP.ItemTokenSize = [26 18];

% Widen the legend after the LaTeX text has been rendered, keeping the right
% edge where it was.
drawnow;
w0 = legP.Position(3);
legP.Position(3) = w0*1.25;
legP.Position(1) = legP.Position(1) - (legP.Position(3) - w0);

exportgraphics(figP, 'co57_sub.pdf', ...
    'ContentType','vector', 'BackgroundColor','white');
fprintf('Exported: co57_sub.pdf  (%.1f x %.1f cm)\n', W_CM, H_CM);
