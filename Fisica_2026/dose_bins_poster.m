%% ======================================================================
%  POSTER VERSION of the two-bin dose method figure (thesis Figure 3)
%
%  Run dose.m first, up to and including the fit section, so that qual, nQ,
%  T_acq, ch2E, k_lo, k_split, k_hi, w and sdw are in the workspace.
%
%  Built at the size it will occupy on the poster, so LaTeX includes it with
%  width=\linewidth and applies almost no scaling. The font sizes below are
%  therefore close to the sizes that get printed.
%
%  Differences from the thesis figure: no grid, thicker rules, larger text,
%  and each weight is printed over the energy range it applies to, in the
%  colour of that bin, so the number and the bin read as one thing.
% ======================================================================

W_CM = 24.0;   H_CM = 24.45;      % final size on the poster, centimetres

FS_TICK = 26;                     % tick labels
FS_LAB  = 30;                     % axis labels
FS_EQ   = 28;                     % the Hp(10) equation, above the axes
FS_BIN  = 40;                     % the N1 and N2 markers inside the fills
FS_VAL  = 26;                     % the weight printed over each bin
FS_QUAL = 24;                     % the N-40 tag in the corner

LW_AXES  = 2.5;
LW_STEP  = 3;                     % histogram outline
LW_SPLIT = 3;                     % bin boundary

col_blue  = [ 26 111 223]/255;    % bin 1
col_green = [ 55 173 107]/255;    % bin 2
col_grey  = [ 81  81  81]/255;

%% ------------------------------------------------- data, N-40 spectrum
q  = nQ;                                  % the hardest quality, N-40
s  = qual(q).ch >= k_lo & qual(q).ch <= min(k_hi, qual(q).k_end);
kk = qual(q).ch(s);
cy = qual(q).cts(s)/T_acq;
el = ch2E(kk-0.5);   eh = ch2E(kk+0.5);

E1 = ch2E(k_lo-0.5);   E2 = ch2E(k_split+0.5);   E3 = ch2E(k_hi+0.5);
i1 = kk <= k_split;    i2 = kk >  k_split;

ymax = max(cy);
ytop = 1.50*ymax;                 % headroom for the two weight labels

%% ------------------------------------------------- figure
figP = figure('Color','w', 'Units','centimeters', ...
              'Position',[2 2 W_CM H_CM], ...
              'PaperUnits','centimeters', 'PaperPosition',[0 0 W_CM H_CM]);
axP = axes(figP);
hold(axP,'on');

fill_bin(axP, el(i1), eh(i1), cy(i1), col_blue,  0.80);
fill_bin(axP, el(i2), eh(i2), cy(i2), col_green, 0.80);

[xv,yv] = step_xy(el, eh, cy);
plot(axP, xv, yv, '-', 'Color', col_grey, 'LineWidth', LW_STEP);

plot(axP, [E2 E2], [0 ytop], '--', 'Color', col_grey, 'LineWidth', LW_SPLIT);

% --- N1 and N2, inside the filled areas -------------------------------
text(axP, mean([E1 E2]), 0.40*ymax, '$N_1$', 'Interpreter','latex', ...
     'FontSize', FS_BIN, 'Color','w', ...
     'HorizontalAlignment','center', 'VerticalAlignment','middle');
text(axP, mean([E2 E3]), 0.40*ymax, '$N_2$', 'Interpreter','latex', ...
     'FontSize', FS_BIN, 'Color','w', ...
     'HorizontalAlignment','center', 'VerticalAlignment','middle');

% --- the weights, each centred over the bin it belongs to --------------
%  Two lines so the label stays inside its own energy range. Bold has to
%  come from the markup, \mathbf and \textbf: the latex interpreter ignores
%  the FontWeight property. If a label ever runs into the other one, drop
%  FS_VAL or shorten the second line.
text(axP, mean([E1 E2]), 1.16*ymax, ...
     {sprintf('$\\mathbf{w_1 = %.0f \\pm %.0f}$', w(1)*1e6, sdw(1)*1e6), ...
      '\textbf{pSv/count}'}, ...
     'Interpreter','latex', 'FontSize', FS_VAL, 'Color', col_blue, ...
     'HorizontalAlignment','center', 'VerticalAlignment','middle');
text(axP, mean([E2 E3]), 1.16*ymax, ...
     {sprintf('$\\mathbf{w_2 = %.0f \\pm %.0f}$', w(2)*1e6, sdw(2)*1e6), ...
      '\textbf{pSv/count}'}, ...
     'Interpreter','latex', 'FontSize', FS_VAL, 'Color', col_green, ...
     'HorizontalAlignment','center', 'VerticalAlignment','middle');

% --- radiation quality, top right corner ------------------------------
text(axP, E3 - 0.02*(E3-E1), 0.96*ytop, '\textbf{N-40}', ...
     'Interpreter','latex', 'FontSize', FS_QUAL, 'Color', col_grey, ...
     'HorizontalAlignment','right', 'VerticalAlignment','top');

hold(axP,'off');

set(axP, ...
    'TickLabelInterpreter','latex', ...
    'FontSize',   FS_TICK, ...
    'LineWidth',  LW_AXES, ...
    'TickDir',    'in', ...
    'TickLength', [0.018 0.010], ...
    'XLim', [E1 E3], 'XTick', 20:5:40, ...
    'YLim', [0 ytop], 'YTick', 0:ceil(ymax/3/50)*50:ytop, ...
    'XMinorTick','off', 'YMinorTick','off', ...
    'Box','on');
grid(axP,'off');                          % no grid on the poster version

xlabel(axP, 'Energy (keV)',           'Interpreter','latex', 'FontSize', FS_LAB);
ylabel(axP, 'Count rate (s$^{-1}$)',  'Interpreter','latex', 'FontSize', FS_LAB);
title(axP,  '$\mathbf{H_p(10) = w_1 N_1 + w_2 N_2}$ \textbf{(Sv)}', ...
      'Interpreter','latex', 'FontSize', FS_EQ, 'Color','k');

set(axP, 'Units','normalized', 'Position', [0.14 0.19 0.83 0.71]);

exportgraphics(figP, 'dose_bins_poster.pdf', ...
    'ContentType','vector', 'BackgroundColor','white');
fprintf('Exported: dose_bins_poster.pdf  (%.1f x %.1f cm)\n', W_CM, H_CM);
fprintf('  bin 1: %.1f to %.1f keV, w1 = %.0f +/- %.0f pSv/count\n', ...
        E1, E2, w(1)*1e6, sdw(1)*1e6);
fprintf('  bin 2: %.1f to %.1f keV, w2 = %.0f +/- %.0f pSv/count\n', ...
        E2, E3, w(2)*1e6, sdw(2)*1e6);

%% ------------------------------------------------- local functions
%  Copies of the helpers in dose.m: local functions are private to the
%  script that defines them, so this file needs its own.
function [xv, yv] = step_xy(el, eh, y)
% Histogram outline from channel edges, so adjacent channels share an edge
% exactly and no gap can open between two filled ranges.
    xv = reshape([el(:)'; eh(:)'], [], 1);
    yv = reshape([y(:)';  y(:)' ], [], 1);
end

function fill_bin(ax, el, eh, y, colour, alpha_face)
    if isempty(el); return; end
    [xv, yv] = step_xy(el, eh, y);
    patch(ax, [xv; xv(end); xv(1)], [yv; 0; 0], colour, ...
          'FaceAlpha', alpha_face, 'EdgeColor','none', 'HandleVisibility','off');
end