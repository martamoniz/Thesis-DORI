%% ======================================================================
%  POSTER VERSION of the energy calibration figure
%
%  Run energy_calibration.m first, up to and including the Co-57 block, so
%  that a, b, se_a, se_b, R2, E_fit, ch_fit, sig_fit, muAm, sigAm, mu_bar,
%  muErr and E_ref are in the workspace.
%
%  The figure is built at the size it will occupy on the poster, so LaTeX
%  includes it with width=\linewidth and applies no scaling. Every font size
%  below is therefore the size that gets printed, and can be compared
%  directly with the poster body text (31 pt).
% ======================================================================

W_CM = 26;   H_CM = 18;          % final size on the poster, centimetres

FS_TICK = 26;                     % tick labels
FS_LAB  = 30;                     % axis labels
FS_LEG  = 22;                     % legend
FS_ANN  = 24;                     % fit-parameter box

LW_AXES = 2.5;                    % axes and box
LW_LINE = 5;                      % fit line
LW_ERR  = 3;                      % error bars
MS_PT   = 14;                     % marker size

blue   = [ 26 111 223]/255;
green  = [ 55 173 107]/255;
yellow = [204 153   0]/255;

E_line  = linspace(18, 62, 200)';
ch_line = a*E_line + b;

figP = figure('Color','w', 'Units','centimeters', ...
              'Position',[2 2 W_CM H_CM], ...
              'PaperUnits','centimeters', 'PaperPosition',[0 0 W_CM H_CM]);
axP = axes(figP);
hold(axP,'on');

plot(axP, E_line, ch_line, '--', 'Color', blue, 'LineWidth', LW_LINE, ...
     'DisplayName', 'Linear Fit');

errorbar(axP, E_fit, ch_fit, sig_fit, 'o', ...
     'Color','k', 'MarkerFaceColor','k', 'MarkerEdgeColor','k', ...
     'MarkerSize', MS_PT, 'LineWidth', LW_ERR, 'CapSize', 14, ...
     'LineStyle','none', 'DisplayName','Bremsstrahlung Endpoints');

errorbar(axP, 59.5, muAm, sigAm, 'o', ...
     'Color', green, 'MarkerFaceColor', green, 'MarkerEdgeColor', green, ...
     'MarkerSize', MS_PT+2, 'LineWidth', LW_ERR, 'CapSize', 14, ...
     'LineStyle','none', 'DisplayName','$^{241}$Am Photopeak');

errorbar(axP, E_ref, mu_bar, muErr, 's', ...
     'Color', yellow, 'MarkerFaceColor', yellow, 'MarkerEdgeColor', yellow, ...
     'MarkerSize', MS_PT+2, 'LineWidth', LW_ERR, 'CapSize', 14, ...
     'LineStyle','none', 'DisplayName','$^{57}$Co Compton Edge');

hold(axP,'off');

set(axP, ...
    'TickLabelInterpreter','latex', ...
    'FontSize',   FS_TICK, ...
    'LineWidth',  LW_AXES, ...
    'TickDir',    'in', ...
    'TickLength', [0.018 0.010], ...
    'XLim', [18 62], 'XTick', 20:10:60, ...   % fewer labels than the thesis
    'YLim', [27 57], 'YTick', 30:5:55, ...
    'XMinorTick','on', 'YMinorTick','off', ...   % no minor ticks on Y
    'Box','on');
grid(axP,'off');                              % no grid on the poster version
axP.XAxis.MinorTickValues = 20:5:60;

xlabel(axP, '$E$ (keV)',        'Interpreter','latex', 'FontSize', FS_LAB);
ylabel(axP, 'Channel (a.u.)',   'Interpreter','latex', 'FontSize', FS_LAB);

legP = legend(axP, 'show');
legP.Interpreter   = 'latex';
legP.FontSize      = FS_LEG;
legP.Location      = 'northwest';
legP.Box           = 'on';
legP.LineWidth     = 1.5;
legP.ItemTokenSize = [40 30];

% MATLAB sizes the legend before the LaTeX interpreter renders the text, so the
% box comes out too narrow for the longest entry. Widen it after a drawnow.
drawnow;
legP.Position(3) = legP.Position(3) * 1.30;
legP.Position(4) = legP.Position(4) * 1.10;

set(axP, 'Units','normalized', 'Position', [0.13 0.15 0.84 0.82]);

% Fit-parameter box. Same construction as the thesis figure: an annotation
% textbox anchored in figure coordinates, text left aligned, box sized to the
% text. Move it by changing the first two numbers of the position vector
% (left edge, bottom edge, both as a fraction of the figure).
annot_str = {'Function: $y = ax + b$', ...
    sprintf('$a = (%.0f \\pm %.0f) \\times 10^{-2}$', a*100, se_a*100), ...
    sprintf('$b = %.0f \\pm %.0f$', b, se_b), ...
    sprintf('$R^2 = %.3f$', R2)};

annotation(figP, 'textbox', [0.52 0.17 0.44 0.30], ...
    'String',              annot_str, ...
    'Interpreter',         'latex', ...
    'FontSize',            FS_ANN, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment',   'middle', ...
    'BackgroundColor',     'white', ...
    'EdgeColor',           blue, ...
    'LineWidth',           2, ...
    'Margin',              8, ...
    'FitBoxToText',        'on');

exportgraphics(figP, 'energy_poster.pdf', ...
    'ContentType','vector', 'BackgroundColor','white');
fprintf('Exported: energy_poster.pdf  (%.1f x %.1f cm)\n', W_CM, H_CM);
