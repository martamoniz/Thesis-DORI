%% ======================================================================
%  POSTER VERSION of the bremsstrahlung spectra with the Duane-Hunt fits
%
%  Top picture of the three-figure panel in RESULTS. Load
%  brem_fit_vs_current.mat first, so that bremData, currents and voltages
%  are in the workspace.
%
%  Built at the size it occupies on the poster, so LaTeX includes it at
%  width=\linewidth and applies no scaling. Every font size below is
%  therefore the size that gets printed.
% ======================================================================

W_CM = 24.0;   H_CM = 12.0;       % final size on the poster, centimetres

FS_TICK = 24;                     % tick labels
FS_LAB  = 27;                     % axis labels
FS_LEG  = 18;                     % legend
FS_TTL  = 26;                     % title

LW_AXES = 2.5;                    % axes and box
LW_LINE = 3.5;                    % spectra
LW_FIT  = 4.5;                    % Duane-Hunt fits

grey = [81 81 81]/255;

I_SEL = 0.250;                    % which tube current to show, mA

%% ------------------------------------------------------------- data
dataPath = 'C:\Users\marta\OneDrive - Universidade de Aveiro\DORI_25-26\RESULTS\18 junho';
load(fullfile(dataPath,'brem_fit_vs_current.mat'));

idxSel = find(abs(currents - I_SEL) < 1e-9, 1);

%% ------------------------------------------------------------- figure
figP = figure('Color','w', 'Units','centimeters', ...
              'Position',[2 2 W_CM H_CM], ...
              'PaperUnits','centimeters', 'PaperPosition',[0 0 W_CM H_CM]);
axP = axes(figP);
hold(axP,'on');

leg_handles = gobjects(0);
leg_labels  = {};

for i = 1:numel(voltages)

    bd = bremData(i, idxSel);

    binsPlot   = bd.bins(bd.startIdx:bd.cutoffIdx);
    countsPlot = movmean(double(bd.corrected(bd.startIdx:bd.cutoffIdx)), 3);

    h = plot(axP, binsPlot, countsPlot, '-', ...
             'Color', bd.color, 'LineWidth', LW_LINE);
    leg_handles(end+1) = h;                                       %#ok<SAGROW>
    leg_labels{end+1}  = sprintf('%d kVp', bd.voltage);           %#ok<SAGROW>

    xFit = linspace(bd.xSel(1), bd.xEnd, 200);
    yFit = polyval(bd.p, xFit);
    yFit(yFit < 0) = 0;
    plot(axP, xFit, yFit, '--', 'Color', grey, 'LineWidth', LW_FIT);
end

h_fit = plot(axP, nan, nan, '--', 'Color', grey, 'LineWidth', LW_FIT);
leg_handles(end+1) = h_fit;
leg_labels{end+1}  = 'Linear Fit';

hold(axP,'off');

set(axP, ...
    'TickLabelInterpreter','latex', ...
    'FontSize',   FS_TICK, ...
    'LineWidth',  LW_AXES, ...
    'TickDir',    'in', ...
    'TickLength', [0.012 0.008], ...
    'XLim', [19 50], 'XTick', 20:5:50, ...
    'XMinorTick','on', 'YMinorTick','off', ...
    'Box','on');
grid(axP,'off');                  % no grid on the poster version
axP.YAxis.Exponent = 4;

xlabel(axP, 'Channel (a.u.)', 'Interpreter','latex', 'FontSize', FS_LAB);
ylabel(axP, 'Counts (a.u.)',  'Interpreter','latex', 'FontSize', FS_LAB);
title(axP,  sprintf('\\textbf{%.3f mA}', I_SEL), ...
      'Interpreter','latex', 'FontSize', FS_TTL);

% The axes are placed before the legend is created, so the legend anchors
% itself inside the final plot box instead of the default one.
set(axP, 'Units','normalized', 'Position', [0.095 0.215 0.885 0.660]);

legP = legend(axP, leg_handles, leg_labels, ...
              'Interpreter','latex', 'Location','northeast', ...
              'Box','on', 'FontSize', FS_LEG, 'NumColumns', 2);
legP.LineWidth     = 1.5;
legP.ItemTokenSize = [30 22];

% MATLAB sizes the legend before the LaTeX interpreter renders the text, so
% the box comes out too narrow. Widen it after a drawnow, keeping the right
% edge where it was.
drawnow;
w0 = legP.Position(3);
legP.Position(3) = w0*1.12;
legP.Position(1) = legP.Position(1) - (legP.Position(3) - w0);

exportgraphics(figP, 'brem_sub.pdf', ...
    'ContentType','vector', 'BackgroundColor','white');
fprintf('Exported: brem_sub.pdf  (%.1f x %.1f cm)\n', W_CM, H_CM);
