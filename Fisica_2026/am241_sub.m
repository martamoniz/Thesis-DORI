%% ======================================================================
%  POSTER VERSION of the Am-241 photopeak spectrum, on the channel axis
%
%  Bottom right picture of the three-figure panel in RESULTS. The channel
%  axis is deliberate: this panel shows the raw spectra, and the arrow that
%  leaves it points at the fit that turns those channels into energies.
%
%  Nothing is refitted. The Gaussian already stored in fitResults is only
%  evaluated on the channel grid.
%
%  The resolution box is drawn only if a and b (the energy calibration) are
%  in the workspace; without them the panel is still complete.
%
%  Built at the size it occupies on the poster, so LaTeX includes it at
%  width=\subwid and applies no scaling.
% ======================================================================

W_CM = 11.6;   H_CM = 10.8;       % final size on the poster, centimetres

FS_TICK = 20;                     % tick labels
FS_LAB  = 23;                     % axis labels
FS_LEG  = 17;                     % legend
FS_TTL  = 23;                     % title
FS_ANN  = 17;                     % resolution

LW_AXES = 2.0;
LW_FIT  = 4;                      % fit curve
MS_PT   = 5;                      % data markers

col_green = [55 173 107]/255;
dark      = [0 0 0];

X_LO = 22;   X_HI = 90;           % channel window drawn
i    = 1;                         % which acquisition to plot

%% ------------------------------------------------------------- data
srcPath = 'C:\Users\marta\OneDrive - Universidade de Aveiro\DORI_25-26\RESULTS\spectrum\1_junho\Sources\';
load([srcPath 'fitResults_Am241.mat'], 'fitResults');

muAll   = [fitResults.mu];
muAm    = mean(muAll);
sigAm   = max(abs(muAll - muAm));
fwhmAll = [fitResults.fwhm];
W       = mean(fwhmAll);
sigW    = max(abs(fwhmAll - W));

gauss_model = @(pG, x) pG(1).*exp(-(x - pG(2)).^2 ./ (2*pG(3).^2)) + pG(4) + pG(5).*x;

bins      = double(fitResults(i).bins);
corrected = double(fitResults(i).corrected);
corrected(corrected < 1) = 1;                 % same floor as the thesis figure
pHat      = fitResults(i).pHat;

xCurve = linspace(40, 65, 500);               % channel range the fit covers
yCurve = gauss_model(pHat, xCurve);
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

h_fit = plot(axP, xCurve, yCurve, '-', 'Color', col_green, 'LineWidth', LW_FIT);

hold(axP,'off');

set(axP, ...
    'TickLabelInterpreter','latex', ...
    'FontSize',   FS_TICK, ...
    'LineWidth',  LW_AXES, ...
    'TickDir',    'in', ...
    'TickLength', [0.018 0.010], ...
    'XLim', [X_LO X_HI], 'XTick', 30:20:90, ...
    'YLim', [0 ytop], ...
    'XMinorTick','on', 'YMinorTick','off', ...
    'Box','on');
grid(axP,'off');                  % no grid on the poster version

xlabel(axP, 'Channel (a.u.)', 'Interpreter','latex', 'FontSize', FS_LAB);
ylabel(axP, 'Counts (a.u.)',  'Interpreter','latex', 'FontSize', FS_LAB);
title(axP,  '\textbf{$^{241}$Am}', 'Interpreter','latex', ...
      'FontSize', FS_TTL, 'Color', dark);

set(axP, 'Units','normalized', 'Position', [0.215 0.195 0.755 0.680]);

legP = legend(axP, h_fit, 'Photopeak Fit', ...
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

% Centroid, width and resolution, in the same boxed form the thesis figure
% uses. It sits under the legend, in the empty upper right of the axes, and
% is anchored in normalised units so it does not move with the data. The
% resolution line is added only if the calibration is in the workspace.
mu_scaled    = round(muAm * 10);
muErr_scaled = round(sigAm * 10);

ann_str = { sprintf('$\\mu = (%d \\pm %d) \\times 10^{-1}$ (a.u.)', ...
                    mu_scaled, muErr_scaled), ...
            sprintf('$FWHM = %.0f \\pm %.0f$ (a.u.)', W, sigW) };

if exist('a','var') && exist('b','var')
    R_res = W/(muAm - b);
    sigR  = R_res*sqrt( (sigW/W)^2 + (sigAm/(muAm - b))^2 );
    ann_str{end+1} = sprintf('$R = %.0f \\pm %.0f$ \\%%', 100*R_res, 100*sigR);
    fprintf('Am-241: R = %.0f +/- %.0f %%\n', 100*R_res, 100*sigR);
end

text(axP, 0.965, 0.735, ann_str, ...
     'Units','normalized', 'Interpreter','latex', 'FontSize', FS_ANN, ...
     'HorizontalAlignment','right', 'VerticalAlignment','top', ...
     'BackgroundColor', [1 1 1], 'EdgeColor', dark, ...
     'LineWidth', 1.2, 'Margin', 5);

exportgraphics(figP, 'am241_sub.pdf', ...
    'ContentType','vector', 'BackgroundColor','white');
fprintf('Exported: am241_sub.pdf  (%.1f x %.1f cm)\n', W_CM, H_CM);
