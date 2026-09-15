%% ======================================================================
%  POSTER VERSION of the dose calibration weight sets
%
%  Same solve as the thesis script: same spectra, same reference Hp(10),
%  three ways of cutting the validated range into bins.
%
%     A   20-25 | 25-40              two bins, boundary at channel 33
%     B   20-30 | 30-40              two bins, boundary at channel 36
%     C   20-25 | 25-30 | 30-40      three bins, square system
%
%  Only the figure block differs. Each panel is built at the size it occupies
%  on the poster, so LaTeX includes it with no scaling and every font size
%  below is the size that gets printed:
%
%     weights2bin25_sub.pdf   11.6 x 10.0 cm   top left of the panel
%     weights2bin30_sub.pdf   11.6 x 10.0 cm   top right
%     weights3bin_sub.pdf     24.0 x 12.0 cm   across the bottom
%
%  No grid, closed box, thicker rules, and the value above each bar is
%  larger and no longer bold. The three panels cannot share pixels per bin,
%  since the poster fixes their widths, so the bar width is set in
%  centimetres instead and comes out identical in all three.
%
%  The filenames are new, so running this does not overwrite the thesis
%  figures in images/.
% ======================================================================

clear; close all; clc;

%% -------------------------------------------------------------- Config
col_blue   = [ 26 111 223]/255;    % bin 1
col_green  = [ 55 173 107]/255;    % bin 2
col_yellow = [204 153   0]/255;    % bin 3
col_bin    = {col_blue, col_green, col_yellow};

T_acq = 120;                 % acquisition time [s], same for all three runs

a_cal = 0.51;                % stored energy calibration: channel = a*E + b
b_cal = 21.0;
ch2E  = @(k) (k - b_cal)/a_cal;

alpha = 0.5537;              % gain compression, k_original = alpha*k + beta
beta  = 12.939;

k_lo = 31;                   % lowest channel of the validated range
k_hi = 41;                   % highest channel (40 keV endpoint of N-40)

qual(1).name='N-25'; qual(1).kVp=25; qual(1).file='25kV_0.500mA_2min.csv';
qual(1).Ka=0.0603; qual(1).uKa=0.0084; qual(1).hpK=0.574;
qual(2).name='N-30'; qual(2).kVp=30; qual(2).file='30kV_0.500mA_2min.csv';
qual(2).Ka=0.2056; qual(2).uKa=0.0196; qual(2).hpK=0.815;
qual(3).name='N-40'; qual(3).kVp=40; qual(3).file='40kV_0.5mA_2min.csv';
qual(3).Ka=0.2820; qual(3).uKa=0.0210; qual(3).hpK=1.210;
nQ = numel(qual);

u_hpK_rel = 0.02;            % ISO 4037-3, 4.1.2: 2 % (k = 1)

% The three schemes. Each row of .edges is [first_channel last_channel], .nom
% holds the nominal energy edges used for the tick labels. .W and .H are the
% size the panel occupies on the poster, in centimetres.
sc(1).edges = [31 33; 34 41];
sc(1).nom   = [20 25 40];
sc(1).file  = 'weights2bin25_sub.pdf';
sc(1).W     = 11.6;   sc(1).H = 10.0;

sc(2).edges = [31 36; 37 41];
sc(2).nom   = [20 30 40];
sc(2).file  = 'weights2bin30_sub.pdf';
sc(2).W     = 11.6;   sc(2).H = 10.0;

sc(3).edges = [31 33; 34 36; 37 41];
sc(3).nom   = [20 25 30 40];
sc(3).file  = 'weights3bin_sub.pdf';
sc(3).W     = 24.0;   sc(3).H = 12.0;
nS = numel(sc);

%% ------------------------------------- Load and compress the spectra
for q = 1:nQ
    Traw = readtable(qual(q).file);
    k = Traw{:,1}; c = Traw{:,2};
    qual(q).ch    = k;
    qual(q).cts   = compress_gain(k, c, alpha, beta);
    qual(q).k_end = floor(a_cal*qual(q).kVp + b_cal);
end

fprintf('Energy calibration: channel = %.2f*E + %.1f   (%.2f keV/channel)\n\n', ...
        a_cal, b_cal, 1/a_cal);

%% ------------------------------------- Reference Hp(10)
H = zeros(nQ,1); uH = zeros(nQ,1);
for q = 1:nQ
    H(q)  = qual(q).Ka * qual(q).hpK;
    uH(q) = H(q)*hypot(qual(q).uKa/qual(q).Ka, u_hpK_rel);
end
sw = 1./uH;

%% ------------------------------------- Solve every scheme
for s = 1:nS
    E  = sc(s).edges;
    nb = size(E,1);

    A = zeros(nQ,nb);
    for q = 1:nQ
        k = qual(q).ch; c = qual(q).cts;
        for j = 1:nb
            m = k >= E(j,1) & k <= min(E(j,2), qual(q).k_end);
            A(q,j) = sum(c(m))/T_acq;
        end
    end

    w    = lsqnonneg(A.*sw, H.*sw);
    Hhat = A*w;
    chi2 = sum(((H-Hhat)./uH).^2);
    Aw   = A.*sw;
    C    = inv(Aw'*Aw)*max(chi2,1);

    sc(s).nb    = nb;
    sc(s).A     = A;
    sc(s).w     = w;
    sc(s).sdw   = sqrt(diag(C));
    sc(s).resp  = Hhat./H;
    sc(s).chi2  = chi2;
    sc(s).dev   = max(abs(Hhat./H - 1));
    sc(s).wfree = (A.*sw)\(H.*sw);          % unconstrained, to expose clipping
    sc(s).nfill = sum(A > 0, 1);            % qualities contributing to each bin
    sc(s).el    = ch2E(E(:,1) - 0.5);
    sc(s).eh    = ch2E(E(:,2) + 0.5);
end

%% ------------------------------------- Console summary
for s = 1:nS
    fprintf('scheme %d   %s\n', s, sprintf('%d ', sc(s).nom));
    for j = 1:sc(s).nb
        fprintf('  bin %d  channels %2d-%2d  %5.1f to %5.1f keV   w = %6.1f +/- %4.1f pSv/count  (%d qualities)\n', ...
                j, sc(s).edges(j,1), sc(s).edges(j,2), sc(s).el(j), sc(s).eh(j), ...
                sc(s).w(j)*1e6, sc(s).sdw(j)*1e6, sc(s).nfill(j));
    end
    fprintf('  response %s\n', mat2str(round(sc(s).resp,4)'));
    fprintf('  max deviation %.2f %%   chi2 = %.3f\n', 100*sc(s).dev, sc(s).chi2);
    if any(abs(sc(s).wfree - sc(s).w) > 1e-12)
        fprintf('  note: unconstrained solution %s pSv/count was clipped at zero\n', ...
                mat2str(round(sc(s).wfree'*1e6,1)));
    end
    if any(sc(s).nfill < 2)
        fprintf('  note: bin(s) %s fed by a single quality, that weight is not over-determined\n', ...
                mat2str(find(sc(s).nfill < 2)));
    end
    fprintf('\n');
end

%% ============================================================ POSTER STYLE
% One bar width in centimetres for all three panels. The poster fixes the
% panel widths, so the bars can only match if the width is set in absolute
% units and converted to x units panel by panel.
BAR_CM = 3.0;

% Font sizes and rules, one set for the two small panels and one for the wide
% one, matching the other poster figures.
for s = 1:nS
    if sc(s).W < 15
        sc(s).FS_TICK = 20;  sc(s).FS_LAB = 23;  sc(s).FS_IN = 20;
        sc(s).LW_AX   = 2.0; sc(s).LW_ERR = 3.0; sc(s).CAP    = 14;
        sc(s).ML = 2.50;  sc(s).MR = 0.50;  sc(s).MB = 2.00;  sc(s).MT = 0.45;
    else
        sc(s).FS_TICK = 24;  sc(s).FS_LAB = 27;  sc(s).FS_IN = 24;
        sc(s).LW_AX   = 2.5; sc(s).LW_ERR = 3.5; sc(s).CAP    = 18;
        sc(s).ML = 3.00;  sc(s).MR = 0.60;  sc(s).MB = 2.40;  sc(s).MT = 0.50;
    end
end

% Panels with the same number of bins share a y axis, so the two on the top
% row can be read against each other. The three-bin panel sits alone
% underneath and gets its own scale.
for s = 1:nS
    m = 0;
    for g = find([sc.nb] == sc(s).nb)
        m = max(m, max(sc(g).w + sc(g).sdw)*1e6);
    end
    [sc(s).ytop, sc(s).ystep] = nice_top(1.12*m);
end

%% ===================================================== FIGURES 1 to 3
for s = 1:nS

    nb    = sc(s).nb;
    xpos  = 1:nb;
    xl    = [0.4 nb+0.6];
    y     = sc(s).w*1e6;
    u     = sc(s).sdw*1e6;
    ytop  = sc(s).ytop;
    ystep = sc(s).ystep;

    axw = sc(s).W - sc(s).ML - sc(s).MR;      % axes box, centimetres
    axh = sc(s).H - sc(s).MB - sc(s).MT;
    bar_w = BAR_CM/(axw/diff(xl));            % bar width in x units

    % The unit label sits one tick-label height below the axis.
    kev_dy = -1.35*sc(s).FS_TICK*2.54/72/axh;

    fig = figure('Color','w', 'Units','centimeters', ...
                 'Position',[2 2 sc(s).W sc(s).H], ...
                 'PaperUnits','centimeters', ...
                 'PaperPosition',[0 0 sc(s).W sc(s).H]);
    ax = axes(fig);
    hold(ax,'on');

    for j = 1:nb
        bar(ax, xpos(j), y(j), bar_w, ...
            'FaceColor', col_bin{j}, 'EdgeColor','none');
    end

    errorbar(ax, xpos, y, u, 'LineStyle','none', 'Color','k', ...
             'LineWidth', sc(s).LW_ERR, 'CapSize', sc(s).CAP);

    % Value above each bar. Bigger than the thesis version and no longer
    % bold: at this size the weight reads without it.
    for j = 1:nb
        text(ax, xpos(j), y(j) + u(j) + 0.035*ytop, sprintf('$%.1f$', y(j)), ...
             'Interpreter','latex', 'FontSize', sc(s).FS_IN, ...
             'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
             'Color','k');
    end

    hold(ax,'off');

    lab = cell(1,nb);
    lab{1} = sprintf('$[%d,%d]$', sc(s).nom(1), sc(s).nom(2));
    for j = 2:nb
        lab{j} = sprintf('$]%d,%d]$', sc(s).nom(j), sc(s).nom(j+1));
    end

    xlim(ax, xl);   ylim(ax, [0 ytop]);

    set(ax, ...
        'TickLabelInterpreter','latex', ...
        'FontSize',   sc(s).FS_TICK, ...
        'LineWidth',  sc(s).LW_AX, ...
        'TickDir',    'in', ...
        'TickLength', [0.015 0.008], ...
        'XTick', xpos, 'XTickLabel', lab, ...
        'XMinorTick','off', 'YMinorTick','on', ...
        'Box','on');
    grid(ax,'off');                       % no grid on the poster version

    yticks(ax, 0:ystep:ytop);
    ax.YAxis.MinorTickValues = 0:ystep/2:ytop;

    ylabel(ax, '$w_j$ (pSv per count)', ...
           'Interpreter','latex', 'FontSize', sc(s).FS_LAB);

    % Unit at the right end of the tick row, in place of an x label.
    tk = text(ax, 0, 0, 'keV', 'Interpreter','latex', ...
              'FontSize', sc(s).FS_TICK, ...
              'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
              'Clipping','off');
    set(tk, 'Units','normalized', 'Position', [1 kev_dy 0]);

    set(ax, 'Units','normalized', ...
            'Position', [sc(s).ML/sc(s).W, sc(s).MB/sc(s).H, ...
                         axw/sc(s).W,      axh/sc(s).H]);

    % Pin the export bounding box with two white hairlines in opposite
    % corners, then write the vector PDF.
    annotation(fig, 'line', [0.001 0.001], [0.001 0.002], 'Color', 'w');
    annotation(fig, 'line', [0.999 0.999], [0.998 0.999], 'Color', 'w');

    exportgraphics(fig, sc(s).file, ...
        'ContentType','vector', 'BackgroundColor','white');
    fprintf('Exported: %-24s (%.1f x %.1f cm, bars %.1f cm)\n', ...
            sc(s).file, sc(s).W, sc(s).H, BAR_CM);
end

%% ------------------------------------- Table out
rows = [];
for s = 1:nS
    for j = 1:sc(s).nb
        rows = [rows; table(s, j, sc(s).edges(j,1), sc(s).edges(j,2), ...
                sc(s).el(j), sc(s).eh(j), sc(s).w(j)*1e6, sc(s).sdw(j)*1e6, ...
                100*sc(s).dev, sc(s).chi2, ...
                'VariableNames',{'scheme','bin','k_first','k_last', ...
                'E_low_keV','E_high_keV','w_pSv_per_count', ...
                'u_w_pSv_per_count','max_dev_pct','chi2'})];  %#ok<AGROW>
    end
end
writetable(rows, 'dose_calibration_binning_schemes.xlsx');

%% ------------------------------------- local functions
function [ytop, step] = nice_top(v)
% Smallest round upper limit that covers v with at most five major ticks.
    cand = [5 10 20 25 50 100 200 250 500 1000];
    for step = cand
        ytop = ceil(v/step)*step;
        if ytop/step <= 5, break; end
    end
end

function out = compress_gain(k, c, alpha, beta)
% Map channel k onto alpha*k + beta and redistribute counts onto the integer
% channel grid. Counts straddling a boundary are split in proportion to the
% overlap, so the total is conserved exactly.
    n   = numel(k);
    out = zeros(n,1);
    lo  = alpha*(k-0.5) + beta;
    hi  = alpha*(k+0.5) + beta;
    for i = 1:n
        if c(i) == 0; continue; end
        j0 = floor(lo(i)+0.5); j1 = floor(hi(i)+0.5);
        for j = j0:j1
            ov = min(hi(i), j+0.5) - max(lo(i), j-0.5);
            if ov > 0 && j >= 0 && j <= n-1
                out(j+1) = out(j+1) + c(i)*ov/(hi(i)-lo(i));
            end
        end
    end
end
