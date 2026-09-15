%% ======================================================================
%  POSTER VERSION of the Am-241 photopeak spectrum, on an energy axis
%
%  Run energy_calibration.m first, up to and including the energy resolution
%  block, so that fitResults, muAm, sigAm, W, sigW, a, b and se_b are in the
%  workspace.
%
%  The thesis figure is drawn against channel number. Here the axis is
%  converted with the stored calibration, E = (channel - b)/a, so the
%  photopeak sits at its nominal 59.5 keV and the FWHM is read in keV.
%  Nothing is refitted: the Gaussian is the one already in fitResults, it is
%  only evaluated on the channel grid and then mapped onto energy.
%
%  The whole spectrum is drawn, but the axis carries ticks only at 20, 30,
%  40, 50 and 60 keV, the validated energy range, so the figure does not put
%  a number on energies the calibration does not cover. The box is drawn by a
%  second, empty axes on top, which is how the frame stays closed without
%  mirroring the ticks onto the top and right edges.
%
%  Built at the size it will occupy on the poster, so LaTeX includes it with
%  width=\linewidth and applies no scaling.
% ======================================================================

W_CM = 24.6;   H_CM = 22.6;       % final size on the poster, centimetres

FS_TICK = 26;                     % tick labels
FS_LAB  = 30;                     % axis labels
FS_LEG  = 24;                     % legend
FS_ANN  = 24;                     % the centroid / FWHM box

LW_AXES = 2.5;
LW_FIT  = 5;                      % fit curve
MS_PT   = 9;                      % data markers

col_green = [ 55 173 107]/255;
dark      = [0 0 0];

E_LO = 0;    E_HI = 110;          % energy window drawn, keV
T_LO = 20;   T_HI = 60;           % validated range: the only ticks numbered
i    = 1;                         % which acquisition to plot

%% ------------------------------------------------- data on an energy axis
ch2E = @(k) (k - b)/a;            % stored calibration, channel -> keV

bins      = double(fitResults(i).bins);
corrected = double(fitResults(i).corrected);
pHat      = fitResults(i).pHat;

gauss_model = @(pG, x) pG(1).*exp(-(x - pG(2)).^2 ./ (2*pG(3).^2)) + pG(4) + pG(5).*x;

Eb = ch2E(bins);
m  = Eb >= E_LO & Eb <= E_HI;
Ex = Eb(m);
Ey = corrected(m);
Ey(Ey < 1) = 1;                   % same floor as the thesis figure

chCurve = linspace(40, 65, 500);  % the channel range the Gaussian was fitted on
Ecurve  = ch2E(chCurve);
ycurve  = gauss_model(pHat, chCurve);

%% ------------------------------------------------- numbers in keV
mu_keV   = ch2E(muAm);
sig_mu   = sigAm/a;
fwhm_keV = W/a;
sig_fwhm = sigW/a;
R_res    = W/(muAm - b);          % FWHM/E, identical in channels or in keV
sigR     = R_res*sqrt( (sigW/W)^2 + (sigAm^2 + se_b^2)/(muAm - b)^2 );

fprintf('Am-241 on the energy axis\n');
fprintf('  centroid   = %.1f +/- %.1f keV   (nominal 59.5 keV)\n', mu_keV, sig_mu);
fprintf('  FWHM       = %.1f +/- %.1f keV\n', fwhm_keV, sig_fwhm);
fprintf('  resolution = %.0f +/- %.0f %%\n', 100*R_res, 100*sigR);

%% ------------------------------------------------- figure
ytop = 3*max(ycurve);

figP = figure('Color','w', 'Units','centimeters', ...
              'Position',[2 2 W_CM H_CM], ...
              'PaperUnits','centimeters', 'PaperPosition',[0 0 W_CM H_CM]);
axP = axes(figP);
hold(axP,'on');

plot(axP, Ex, Ey, 'o', 'MarkerSize', MS_PT, ...
     'MarkerEdgeColor', dark, 'MarkerFaceColor', dark, ...
     'HandleVisibility','off');

h_fit = plot(axP, Ecurve, ycurve, '-', 'Color', col_green, 'LineWidth', LW_FIT);

hold(axP,'off');

set(axP, ...
    'TickLabelInterpreter','latex', ...
    'FontSize',   FS_TICK, ...
    'LineWidth',  LW_AXES, ...
    'TickDir',    'in', ...
    'TickLength', [0.018 0.010], ...
    'XLim', [E_LO E_HI], 'XTick', T_LO:10:T_HI, ...
    'YLim', [0 ytop], ...
    'XMinorTick','off', 'YMinorTick','off', ...
    'Box','off');                 % the frame is drawn separately, see below
grid(axP,'off');                  % no grid on the poster version

xlabel(axP, 'Energy (keV)',  'Interpreter','latex', 'FontSize', FS_LAB);
ylabel(axP, 'Counts (a.u.)', 'Interpreter','latex', 'FontSize', FS_LAB);
title(axP,  '\textbf{$^{241}$Am}', 'Interpreter','latex', 'FontSize', 28, 'Color', dark);

legP = legend(axP, h_fit, 'Photopeak Fit', ...
              'Interpreter','latex', 'Location','northeast', ...
              'Box','on', 'FontSize', FS_LEG);
legP.LineWidth     = 1.5;
legP.ItemTokenSize = [40 30];

% MATLAB measures the legend before the LaTeX interpreter has rendered the
% text, so the box comes out too narrow and the last word spills over the
% edge. Widen it after a drawnow, keeping the right edge where it was.
drawnow;
w0 = legP.Position(3);
legP.Position(3) = w0*1.35;
legP.Position(1) = legP.Position(1) - (legP.Position(3) - w0);
legP.Position(4) = legP.Position(4)*1.10;

set(axP, 'Units','normalized', 'Position', [0.15 0.19 0.82 0.72]);

% Closed frame with no ticks of its own: an empty axes sitting exactly on top
% of the plot. Ticks stay on the left and bottom edges only.
axFrame = axes('Parent', figP, 'Units','normalized', 'Position', axP.Position, ...
               'Box','on', 'XTick',[], 'YTick',[], 'Color','none', ...
               'XLim',[0 1], 'YLim',[0 1], 'LineWidth', LW_AXES, ...
               'HandleVisibility','off', 'HitTest','off');

% Centroid and width. Move the box by changing the first two numbers of the
% position vector (left edge, bottom edge, fractions of the figure); the box
% sizes itself to the text.
annot_str = { sprintf('$\\mu = %.1f \\pm %.1f$ keV', mu_keV, sig_mu), ...
              sprintf('FWHM $= %.0f \\pm %.0f$ keV', fwhm_keV, sig_fwhm) };

annotation(figP, 'textbox', [0.38 0.45 0.4 0.28], ...
    'String',              annot_str, ...
    'Interpreter',         'latex', ...
    'FontSize',            FS_ANN, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment',   'middle', ...
    'BackgroundColor',     'white', ...
    'EdgeColor',           col_green, ...
    'LineWidth',           2, ...
    'Margin',              8, ...
    'FitBoxToText',        'on');

exportgraphics(figP, 'am241_poster.pdf', ...
    'ContentType','vector', 'BackgroundColor','white');
fprintf('Exported: am241_poster.pdf  (%.1f x %.1f cm)\n', W_CM, H_CM);