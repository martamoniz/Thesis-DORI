import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

try:
    from scipy.signal import savgol_filter
except Exception:
    savgol_filter = None

# =========================
# CONFIG
# =========================
base_path = r"C:\Users\marta\OneDrive - Universidade de Aveiro\DORI_25-26\wavesurfer 3054\ADC ESP32"

files = {
    "shaper":  "SHAPER.xlsx",
    "trigger": "TRIGGER.xlsx",
    "ph":      "P&H.xlsx",
}

# target horizontal resolution (controls downsampling quality)
N_PIXELS = 1200

# optional smoothing (applied before downsampling) – for shaper only
ENABLE_SMOOTHING = True
SMOOTH_WINDOW    = 21
SMOOTH_POLYORDER = 3
SMOOTH_SHAPER    = True

# axis limits
USE_FIXED_LIMITS = True
SHAPER_YLIM  = (-0.1, 0.6)
TRIGGER_YLIM = (-0.5, 3.8)   # right axis – trigger goes to 3.3 V
PH_YLIM      = (-0.1, 0.6)   # left axis – same scale as shaper

# ── Step-aware downsampling parameters ───────────────────────────────────────
# A sample is "the same plateau" as its neighbour if the absolute difference
# is below this threshold.  Tune per signal if needed.
TRIGGER_FLAT_TOL = 0.05   # V  (very tight – trigger is nearly ideal square wave)
PH_FLAT_TOL      = 0.02   # V  (P&H has slow droop, keep a bit of tolerance)

# =========================
# LOAD FUNCTION
# =========================
def load_xlsx(path):
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Excel file not found: {path}\n"
            "Check base_path and filenames."
        )
    try:
        data = pd.read_excel(path)
    except PermissionError as exc:
        raise PermissionError(
            f"Permission denied reading: {path}\n"
            "Close the file if it is open, or move it out of the OneDrive sync folder."
        ) from exc

    t = data.iloc[:, 0].values.astype(float)
    y = data.iloc[:, 1].values.astype(float)
    return t, y


# =========================
# LARGEST-TRIANGLE-THREE-BUCKETS (LTTB) – for smooth analogue signals
# =========================
def downsample_lttb(t, y, n_out):
    n = len(t)
    if n_out >= n or n_out < 3:
        return t, y

    bucket_size = (n - 2) / (n_out - 2)
    a_idx = 0
    t_out = [t[0]]
    y_out = [y[0]]

    for i in range(n_out - 2):
        start    = int(np.floor((i + 1) * bucket_size)) + 1
        end      = min(int(np.floor((i + 2) * bucket_size)) + 1, n)
        avg_start = int(np.floor((i + 2) * bucket_size)) + 1
        avg_end   = min(int(np.floor((i + 3) * bucket_size)) + 1, n)

        avg_x = np.mean(t[avg_start:avg_end]) if avg_start < avg_end else t[-1]
        avg_y = np.mean(y[avg_start:avg_end]) if avg_start < avg_end else y[-1]

        seg_t = t[start:end]
        seg_y = y[start:end]
        if seg_t.size == 0:
            continue

        area  = np.abs((t[a_idx] - avg_x) * (seg_y - y[a_idx])
                       - (t[a_idx] - seg_t) * (avg_y - y[a_idx]))
        idx   = int(np.argmax(area))
        a_idx = start + idx
        t_out.append(t[a_idx])
        y_out.append(y[a_idx])

    t_out.append(t[-1])
    y_out.append(y[-1])
    return np.array(t_out), np.array(y_out)


# =========================
# STEP-AWARE DOWNSAMPLING  – for digital / plateau signals
# =========================
def downsample_step_aware(t, y, flat_tol, n_out=None):
    """
    Preserves every transition edge exactly.
    Within flat plateaux, keeps only the first and last sample of each run
    (so the plot draws a perfect horizontal line there).
    If the result still has more than n_out points it is further thinned by
    keeping every k-th interior point of each plateau – transitions are NEVER
    removed.

    Parameters
    ----------
    t, y      : raw arrays
    flat_tol  : |Δy| below which consecutive samples are considered "same level"
    n_out     : optional cap on output points (None = no cap beyond edge compression)
    """
    if len(t) < 3:
        return t, y

    # ── 1.  Find transition indices ──────────────────────────────────────────
    diff       = np.abs(np.diff(y))
    is_edge    = diff > flat_tol          # True where a step happens
    # An edge between i and i+1 → keep both i and i+1
    edge_mask  = np.zeros(len(t), dtype=bool)
    edge_mask[:-1] |= is_edge
    edge_mask[1:]  |= is_edge

    # Always keep first and last
    edge_mask[0]  = True
    edge_mask[-1] = True

    # ── 2.  Within each plateau, keep first + last sample of the run ─────────
    #  Walk the array and mark the first/last of every contiguous flat run.
    plateau_mask = np.zeros(len(t), dtype=bool)
    in_plateau   = False
    run_start    = 0

    for i in range(len(t)):
        if edge_mask[i]:
            if in_plateau and run_start != i:
                plateau_mask[run_start] = True
                plateau_mask[i - 1]     = True
            in_plateau = False
        else:
            if not in_plateau:
                run_start  = i
                in_plateau = True

    if in_plateau:
        plateau_mask[run_start] = True
        plateau_mask[-1]        = True

    keep = edge_mask | plateau_mask

    t_k = t[keep]
    y_k = y[keep]

    # ── 3.  Optional cap – thin interior plateau points if still too many ────
    if n_out is not None and len(t_k) > n_out:
        # Identify which kept points are edges vs interior plateau
        diff2    = np.abs(np.diff(y_k))
        is_edge2 = np.zeros(len(t_k), dtype=bool)
        is_edge2[:-1] |= diff2 > flat_tol
        is_edge2[1:]  |= diff2 > flat_tol
        is_edge2[0]    = True
        is_edge2[-1]   = True

        interior = np.where(~is_edge2)[0]
        edges    = np.where( is_edge2)[0]

        # How many interior points can we keep?
        budget  = max(0, n_out - len(edges))
        if budget < len(interior):
            keep2         = np.zeros(len(t_k), dtype=bool)
            keep2[edges]  = True
            chosen        = interior[np.round(
                np.linspace(0, len(interior) - 1, budget)
            ).astype(int)]
            keep2[chosen] = True
            t_k = t_k[keep2]
            y_k = y_k[keep2]

    return t_k, y_k


# =========================
# TRIGGER-SPECIFIC DOWNSAMPLING
# =========================
def downsample_trigger(t, y, flat_tol=0.05):
    """
    Special case for a near-perfect square wave:
      • Compute global median of LOW samples and HIGH samples
      • Replace every plateau segment with exactly its median level
      • Return only the transition edges + one point per plateau
    This gives perfectly horizontal lines with zero noise.
    """
    # ── Detect high / low by kmeans-style split ──────────────────────────────
    midpoint  = (np.max(y) + np.min(y)) / 2
    low_mask  = y <= midpoint
    high_mask = ~low_mask

    low_level  = np.median(y[low_mask])  if np.any(low_mask)  else np.min(y)
    high_level = np.median(y[high_mask]) if np.any(high_mask) else np.max(y)

    # ── Snap every sample to its level ──────────────────────────────────────
    y_clean = np.where(high_mask, high_level, low_level)

    # ── Keep only transition points (level changes) ──────────────────────────
    # Include the sample before AND after every edge, plus first/last
    changes   = np.where(np.diff(y_clean) != 0)[0]
    keep      = np.zeros(len(t), dtype=bool)
    keep[0]   = True
    keep[-1]  = True
    keep[changes]     = True
    keep[changes + 1] = True

    return t[keep], y_clean[keep]


# =========================
# SMOOTH FUNCTION (shaper only)
# =========================
def smooth_signal(y, enable, window, polyorder):
    if not enable:
        return y
    if window < 5:
        return y
    if window % 2 == 0:
        window += 1
    if window >= len(y):
        window = len(y) if len(y) % 2 == 1 else len(y) - 1
    if window < 5 or window >= len(y):
        return y
    if polyorder >= window:
        polyorder = max(2, window - 2)

    if savgol_filter is not None:
        return savgol_filter(y, window_length=window, polyorder=polyorder)

    kernel = np.ones(window, dtype=float) / window
    return np.convolve(y, kernel, mode="same")


# =========================
# LOAD DATA
# =========================
t_s,  y_s  = load_xlsx(os.path.join(base_path, files["shaper"]))
t_tr, y_tr = load_xlsx(os.path.join(base_path, files["trigger"]))
t_ph, y_ph = load_xlsx(os.path.join(base_path, files["ph"]))

# =========================
# PROCESS SHAPER  (smooth → LTTB)
# =========================
y_s_sm    = smooth_signal(y_s, ENABLE_SMOOTHING and SMOOTH_SHAPER,
                          SMOOTH_WINDOW, SMOOTH_POLYORDER)
t_s_d, y_s_d = downsample_lttb(t_s, y_s_sm, N_PIXELS)

# =========================
# PROCESS TRIGGER  (snap to levels → keep edges only)
# =========================
t_tr_d, y_tr_d = downsample_trigger(t_tr, y_tr, flat_tol=TRIGGER_FLAT_TOL)

# =========================
# PROCESS P&H  (step-aware compression)
# =========================
t_ph_d, y_ph_d = downsample_step_aware(t_ph, y_ph,
                                        flat_tol=PH_FLAT_TOL,
                                        n_out=N_PIXELS)

# =========================
# PLOT
# =========================
fig, ax_left = plt.subplots(figsize=(9, 4.5))

# Right axis for Trigger (3.3 V scale)
ax_right = ax_left.twinx()

# ── Shaper (left axis) ───────────────────────────────────────────────────────
l1, = ax_left.plot(t_s_d,  y_s_d,  linewidth=0.8, color="C0", label="Shaper")

# ── P&H (left axis) ──────────────────────────────────────────────────────────
l2, = ax_left.plot(t_ph_d, y_ph_d, linewidth=0.8, color="C2", label="P&H")

# ── Trigger (right axis) ─────────────────────────────────────────────────────
l3, = ax_right.plot(t_tr_d, y_tr_d, linewidth=1.0, color="C1",
                    linestyle="--", label="Trigger", alpha=0.85)

# ── Axis limits ──────────────────────────────────────────────────────────────
if USE_FIXED_LIMITS:
    ax_left.set_ylim(*SHAPER_YLIM)    # covers both Shaper and P&H
    ax_right.set_ylim(*TRIGGER_YLIM)

# ── Labels ───────────────────────────────────────────────────────────────────
ax_left.set_xlabel("Time (s)")
ax_left.set_ylabel("Amplitude (V)  —  Shaper / P&H")
ax_right.set_ylabel("Amplitude (V)  —  Trigger", color="C1")
ax_right.tick_params(axis="y", labelcolor="C1")

# ── Grid (from left axis only, to avoid double grid) ────────────────────────
ax_left.grid(True, linestyle=":", linewidth=0.5, alpha=0.7)

# ── Combined legend ───────────────────────────────────────────────────────────
ax_left.legend(handles=[l1, l2, l3], loc="upper right")

# =========================
# EXPORT
# =========================
plt.tight_layout()
plt.savefig("waveforms_publication.pdf", dpi=300)
plt.savefig("waveforms_publication.png", dpi=300)
plt.show()
