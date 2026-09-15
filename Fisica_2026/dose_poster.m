%% ======================================================================
%  POSTER VERSION of the dose calibration closure figure (Figure 7)
%
%  Run dose.m first, up to and including the fit section, so that resp, loo,
%  qual and nQ are in the workspace.
%
%  Built at the size it will occupy on the poster, so LaTeX includes it with
%  width=\linewidth and applies no scaling. The font sizes below are the
%  sizes that get printed, comparable with the poster body text (31 pt).
% ======================================================================

W_CM = 26;   H_CM = 15;          % final size on the poster, centimetres

FS_TICK = 26;                     % tick labels
FS_LAB  = 30;                     % axis labels
FS_LEG  = 24;                     % legend
FS_VAL  = 24;                     % value labels on the bars

LW_AXES = 2.5;
LW_ZERO = 2.5;                    % zero line

col_blue  = [ 26 111 223]/255;
col_green = [ 55 173 107]/255;
col_grey  = [ 81  81  81]/255;

d_fit = 100*(resp(:) - 1);        % signed fit residual, per cent
d_loo = 100*loo(:);               % signed leave-one-out error, per cent

ybot = -0.45;   ytop = 0.35;      % a little headroom for the value labels

figP = figure('Color','w', 'Units','centimeters', ...
              'Position',[2 2 W_CM H_CM], ...
              'PaperUnits','centimeters', 'PaperPosition',[0 0 W_CM H_CM]);
axP = axes(figP);
hold(axP,'on');

plot(axP, [0.4 nQ+0.6], [0 0], '-', 'Color', col_grey, ...
     'LineWidth', LW_ZERO, 'HandleVisibility','off');

hb = bar(axP, 1:nQ, [d_fit d_loo], 1, 'EdgeColor','none');
hb(1).FaceColor = col_blue;
hb(2).FaceColor = col_green;

off = 0.022*(ytop - ybot);        % label gap above or below the bar tip
for j = 1:2
    v  = hb(j).YData;  xe = hb(j).XEndPoints;
    for q = 1:nQ
        if v(q) >= 0, s = +1; va = 'bottom'; else, s = -1; va = 'top'; end
        text(axP, xe(q), v(q) + s*off, sprintf('$%+.2f$', v(q)), ...
             'Interpreter','latex', 'FontSize', FS_VAL, ...
             'HorizontalAlignment','center', 'VerticalAlignment', va, ...
             'Color','k');
    end
end
hold(axP,'off');

set(axP, ...
    'TickLabelInterpreter','latex', ...
    'FontSize',   FS_TICK, ...
    'LineWidth',  LW_AXES, ...
    'TickDir',    'in', ...
    'TickLength', [0.018 0.010], ...
    'XLim', [0.4 nQ+0.6], 'XTick', 1:nQ, 'XTickLabel', {qual.name}, ...
    'XMinorTick','off', ...
    'YLim', [ybot ytop], 'YTick', -0.4:0.2:0.2, ...
    'YMinorTick','on', ...
    'Box','on');
grid(axP,'off');                              % no grid on the poster version
axP.YAxis.MinorTickValues = -0.4:0.1:0.3;

xlabel(axP, 'Radiation Quality', 'Interpreter','latex', 'FontSize', FS_LAB);
ylabel(axP, 'Deviation (\%)',    'Interpreter','latex', 'FontSize', FS_LAB);

legP = legend(axP, hb, {'Fit Residual','Leave-One-Out'}, ...
              'Interpreter','latex', 'Location','northeast', ...
              'Box','on', 'FontSize', FS_LEG);
legP.LineWidth     = 1.5;
legP.ItemTokenSize = [40 30];

set(axP, 'Units','normalized', 'Position', [0.14 0.20 0.83 0.77]);

exportgraphics(figP, 'dose_calibration_poster.pdf', ...
    'ContentType','vector', 'BackgroundColor','white');
fprintf('Exported: dose_calibration_poster.pdf  (%.1f x %.1f cm)\n', W_CM, H_CM);
