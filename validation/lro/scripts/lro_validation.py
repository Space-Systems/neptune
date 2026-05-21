"""
LRO SPICE validation workflow for NEPTUNE Moon-centered propagation.

Workflow
--------
  1. Extract LRO initial state from SPICE at --start
  2. Patch neptune.inp  (input type 6, epochs, state vector, enable MCRF outputs)
  3. Optionally run NEPTUNE  (--run)
  4. Compare neptune.mcs / neptune.osm against SPICE truth and report residuals

Required SPICE kernels (set paths in CONFIG):
  naif0012.tls     leap-second kernel
  de421.bsp        planetary ephemeris  (same DE NEPTUNE uses internally)
  lrorg_*.bsp      LRO reconstructed SPK  (NAIF PDS)

Usage examples
--------------
  # Extract state at 2009-07-01, propagate 7 days, run NEPTUNE, compare:
  python lro_validation.py --start "2009-07-01T00:00:00" --duration 7 --run

  # Use already-run output, skip patching / running:
  python lro_validation.py --no-patch --plot

  # Save plots to PNG:
  python lro_validation.py --start "2009-07-01T00:00:00" --run --save-plots
"""

import argparse
import re
import subprocess
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import spiceypy as spice
except ImportError:
    sys.exit("spiceypy not found.  Install with:  pip install spiceypy")

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

# ---------------------------------------------------------------------------
# CONFIG — adjust to your local layout
# ---------------------------------------------------------------------------
WORK_DIR    = Path(__file__).parent.parent / "work"
KERNEL_DIR  = WORK_DIR / "data"
INP_FILE    = WORK_DIR / "input"  / "neptune.inp"
OUT_DIR     = WORK_DIR / "output"
NEPTUNE_CMD = [str(WORK_DIR / "start_neptune.sh")]  # launched from WORK_DIR

# https://naif.jpl.nasa.gov/pub/naif/pds/data/lro-l-spice-6-v1.0/lrosp_1000/aareadme.htm
KERNELS = [
    KERNEL_DIR / "naif0012.tls",
    KERNEL_DIR / "de421.bsp",
    # Add LRO SPK files that cover your window, e.g.:
    KERNEL_DIR / "lrorg_2025258_2025349_v01.bsp",
]

LRO_NAIF_ID  = "-85"   # LRO
MOON_NAIF_ID = "301"   # Moon
SPICE_FRAME  = "J2000" # = GCRF-aligned; matches NEPTUNE MCRF axes exactly

MOON_MU = 4902.800066  # km^3/s^2 — DE-430 (mirrors slam_moon_astro.f90)
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Moon Cartesian -> Keplerian  (Python mirror of moon_rv2coe in libslam)
# ---------------------------------------------------------------------------
def moon_rv2coe(r, v):
    """
    Selenocentric osculating Keplerian elements from Cartesian state.

    Parameters
    ----------
    r : (3,)  position  [km]
    v : (3,)  velocity  [km/s]

    Returns
    -------
    dict  sma [km], ecc [-], inc/raan/aop/tran/man [deg]
    """
    r, v = np.asarray(r, float), np.asarray(v, float)
    mu   = MOON_MU
    eps  = 1e-9

    r_mag = np.linalg.norm(r)
    v_mag = np.linalg.norm(v)
    h     = np.cross(r, v)
    h_mag = np.linalg.norm(h)

    if h_mag < eps:
        return {k: np.nan for k in ("sma", "ecc", "inc", "raan", "aop", "tran", "man")}

    n     = np.array([-h[1], h[0], 0.0])
    n_mag = np.linalg.norm(n)
    rdotv = float(np.dot(r, v))
    e_vec = ((v_mag**2 - mu / r_mag) * r - rdotv * v) / mu
    ecc   = np.linalg.norm(e_vec)

    sme = 0.5 * v_mag**2 - mu / r_mag
    sma = (-mu / (2.0 * sme)) if abs(sme) > eps else np.inf
    inc = np.degrees(np.arccos(np.clip(h[2] / h_mag, -1, 1)))

    raan = (np.degrees(np.arccos(np.clip(n[0] / n_mag, -1, 1)))
            if n_mag > eps else np.nan)
    if n_mag > eps and n[1] < 0:
        raan = 360.0 - raan

    if ecc > eps and n_mag > eps:
        aop = np.degrees(np.arccos(np.clip(np.dot(n, e_vec) / (n_mag * ecc), -1, 1)))
        if e_vec[2] < 0:
            aop = 360.0 - aop
    else:
        aop = np.nan

    if ecc > eps:
        tran = np.degrees(np.arccos(np.clip(np.dot(e_vec, r) / (ecc * r_mag), -1, 1)))
        if rdotv < 0:
            tran = 360.0 - tran
    else:
        tran = np.nan

    if not np.isnan(tran) and ecc < 1.0:
        nu   = np.radians(tran)
        sinE = np.sqrt(1 - ecc**2) * np.sin(nu) / (1 + ecc * np.cos(nu))
        cosE = (ecc + np.cos(nu)) / (1 + ecc * np.cos(nu))
        E    = np.arctan2(sinE, cosE)
        man  = np.degrees(E - ecc * np.sin(E)) % 360.0
    else:
        man = np.nan

    return dict(sma=sma, ecc=ecc, inc=inc, raan=raan, aop=aop, tran=tran, man=man)


# ---------------------------------------------------------------------------
# SPICE helpers
# ---------------------------------------------------------------------------
def load_kernels():
    missing = [str(k) for k in KERNELS if not Path(k).exists()]
    if missing:
        sys.exit(
            "Missing SPICE kernels:\n  " + "\n  ".join(missing) +
            f"\nPlace them in {KERNEL_DIR} or update KERNELS in this script."
        )
    for k in KERNELS:
        spice.furnsh(str(k))
    print(f"Loaded {len(KERNELS)} SPICE kernel(s).")


def iso_to_et(iso_str):
    """Convert ISO-8601 UTC string to SPICE ET."""
    return spice.utc2et(iso_str.replace("T", " ") + " UTC")


def et_to_dt(et):
    """Convert SPICE ET to Python datetime (UTC)."""
    utc = spice.et2utc(et, "ISOC", 3)   # "YYYY-MM-DDTHH:MM:SS.mmm"
    return datetime.strptime(utc, "%Y-%m-%dT%H:%M:%S.%f")


def neptune_date_to_et(date_str):
    """Convert NEPTUNE output date string (YYYY-MM-DDThh:mm:ss.xxxxxxZ) to ET."""
    utc = date_str.rstrip("Z").replace("T", " ") + " UTC"
    return spice.utc2et(utc)


def lro_state_moon_j2000(et):
    """LRO state [r km, v km/s] w.r.t. Moon centre in J2000 (= MCRF axes)."""
    state, _ = spice.spkezr(LRO_NAIF_ID, et, SPICE_FRAME, "NONE", MOON_NAIF_ID)
    return np.asarray(state)


# ---------------------------------------------------------------------------
# Step 1 — extract LRO initial state
# ---------------------------------------------------------------------------
def extract_initial_state(start_iso):
    """
    Return (r_mcrf [km], v_mcrf [km/s], et, dt) for LRO at start_iso (UTC).
    """
    et = iso_to_et(start_iso)
    state = lro_state_moon_j2000(et)
    r, v  = state[:3], state[3:]
    dt    = et_to_dt(et)

    kep = moon_rv2coe(r, v)
    print(f"\nLRO initial state at {dt.isoformat()} UTC  (MCRF / J2000, Moon centre)")
    print(f"  r  = [{r[0]:14.6f}, {r[1]:14.6f}, {r[2]:14.6f}]  km")
    print(f"  v  = [{v[0]:12.9f}, {v[1]:12.9f}, {v[2]:12.9f}]  km/s")
    print(f"  sma={kep['sma']:.3f} km   ecc={kep['ecc']:.6f}   inc={kep['inc']:.4f} deg")

    return r, v, et, dt


# ---------------------------------------------------------------------------
# Step 2 — patch neptune.inp
# ---------------------------------------------------------------------------
def _first_value_lines_after(lines, anchor, count):
    """
    Return line indices of the first `count` non-comment, non-blank lines
    that follow the first line containing `anchor` as a comment.
    """
    found_anchor = False
    result = []
    for i, raw in enumerate(lines):
        s = raw.strip()
        if not found_anchor:
            if s.startswith("#") and anchor.lower() in s.lower():
                found_anchor = True
        else:
            if s and not s.startswith("#"):
                result.append(i)
                if len(result) >= count:
                    break
    return result


def patch_neptune_inp(inp_path, start_dt, end_dt, r_mcrf, v_mcrf):
    """
    Update neptune.inp in-place:
      - input type  → 6  (MCRF Cartesian)
      - begin/end epochs
      - Cartesian state vector (MCRF components)
      - enable OUTPUT_CSV_MCRF and OUTPUT_OSC_MCRF
    """
    lines = Path(inp_path).read_text().splitlines(keepends=True)

    # Helper: replace the content of a specific line index
    def set_line(idx, new_content):
        lines[idx] = new_content if new_content.endswith("\n") else new_content + "\n"

    # --- input type --------------------------------------------------------
    idxs = _first_value_lines_after(lines, "Input type (state vector)", 1)
    if not idxs:
        print("WARNING: could not find 'Input type (state vector)' anchor in neptune.inp")
    else:
        set_line(idxs[0], "  6\n")

    # --- epochs ------------------------------------------------------------
    def fmt_epoch(dt):
        return (f"  {dt.year:04d} {dt.month:02d} {dt.day:02d} "
                f"{dt.hour:02d} {dt.minute:02d} {dt.second:02d}"
                f".{dt.microsecond // 1000:03d}\n")

    idxs = _first_value_lines_after(lines, "Propagation time span", 2)
    if len(idxs) < 2:
        print("WARNING: could not find epoch lines in neptune.inp")
    else:
        set_line(idxs[0], fmt_epoch(start_dt))
        set_line(idxs[1], fmt_epoch(end_dt))

    # --- Cartesian state vector --------------------------------------------
    idxs = _first_value_lines_after(lines, "Cartesian state vector (if input type", 6)
    if len(idxs) < 6:
        print("WARNING: could not find all 6 state-vector lines in neptune.inp")
    else:
        comps = list(r_mcrf) + list(v_mcrf)
        labels = [
            "radius   x-component (km)",
            "radius   y-component (km)",
            "radius   z-component (km)",
            "velocity x-component (km/s)",
            "velocity y-component (km/s)",
            "velocity z-component (km/s)",
        ]
        for i, (val, lbl) in enumerate(zip(comps, labels)):
            set_line(idxs[i], f"  {val:.9f}               {lbl}\n")

    # --- enable MCRF output switches --------------------------------------
    for i, raw in enumerate(lines):
        if "(MCRF, Moon-centred)" in raw and not raw.strip().startswith("#"):
            lines[i] = re.sub(r"^\s*[01]", "  1", raw)

    Path(inp_path).write_text("".join(lines))
    print(f"\nPatched {inp_path}")
    print(f"  Input type : 6  (MCRF Cartesian)")
    print(f"  Start      : {start_dt.isoformat()}")
    print(f"  End        : {end_dt.isoformat()}")
    print(f"  r_mcrf     : {r_mcrf}")
    print(f"  v_mcrf     : {v_mcrf}")
    print(f"  MCRF outputs enabled.")


# ---------------------------------------------------------------------------
# Step 3 — run NEPTUNE  (optional)
# ---------------------------------------------------------------------------
def run_neptune():
    print(f"\nRunning NEPTUNE from {WORK_DIR} ...")
    try:
        subprocess.run(
            NEPTUNE_CMD,
            cwd=str(WORK_DIR),
            capture_output=False,
            check=True,
        )
        print("NEPTUNE finished successfully.")
    except FileNotFoundError:
        sys.exit(f"NEPTUNE command not found: {NEPTUNE_CMD}\nUpdate NEPTUNE_CMD in this script.")
    except subprocess.CalledProcessError as exc:
        sys.exit(f"NEPTUNE exited with code {exc.returncode}.")


# ---------------------------------------------------------------------------
# Step 4 — read NEPTUNE output files
# ---------------------------------------------------------------------------
def read_neptune_mcs(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            tok = line.split()
            if len(tok) < 8:
                continue
            rows.append(dict(
                date=tok[0], mjd=float(tok[1]),
                rx=float(tok[2]), ry=float(tok[3]), rz=float(tok[4]),
                vx=float(tok[5]), vy=float(tok[6]), vz=float(tok[7]),
            ))
    return pd.DataFrame(rows)


def read_neptune_osm(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            tok = line.split()
            if len(tok) < 9:
                continue
            rows.append(dict(
                date=tok[0], mjd=float(tok[1]),
                sma=float(tok[2]),  ecc=float(tok[3]),
                inc=float(tok[4]),  raan=float(tok[5]),
                aop=float(tok[6]),  tran=float(tok[7]), man=float(tok[8]),
            ))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Step 4 — compare Cartesian
# ---------------------------------------------------------------------------
def compare_mcs(mcs_path):
    df = read_neptune_mcs(mcs_path)
    print(f"\nRead {len(df)} epochs from {mcs_path}")

    records = []
    for _, row in df.iterrows():
        try:
            et = neptune_date_to_et(row["date"])
            sp = lro_state_moon_j2000(et)
        except spice.utils.exceptions.SpiceyError as exc:
            print(f"  SPICE error at {row['date']}: {exc}")
            continue

        nep_r = np.array([row["rx"], row["ry"], row["rz"]])
        nep_v = np.array([row["vx"], row["vy"], row["vz"]])
        dr    = nep_r - sp[:3]
        dv    = nep_v - sp[3:]

        records.append(dict(
            mjd=row["mjd"],
            dr_x=dr[0], dr_y=dr[1], dr_z=dr[2],
            dv_x=dv[0], dv_y=dv[1], dv_z=dv[2],
            dr_m=np.linalg.norm(dr) * 1e3,     # m
            dv_mps=np.linalg.norm(dv) * 1e3,   # mm/s
        ))

    res = pd.DataFrame(records)
    print("\n--- Cartesian residuals  (NEPTUNE - SPICE/LRO) ---")
    print(f"  Position RMS : {_rms(res['dr_m']):.3f} m")
    print(f"  Position max : {res['dr_m'].max():.3f} m")
    print(f"  Velocity RMS : {_rms(res['dv_mps']):.3f} mm/s")
    print(f"  Velocity max : {res['dv_mps'].max():.3f} mm/s")
    return res


# ---------------------------------------------------------------------------
# Step 4 — compare Keplerian
# ---------------------------------------------------------------------------
def compare_osm(osm_path):
    df = read_neptune_osm(osm_path)
    print(f"\nRead {len(df)} epochs from {osm_path}")

    records = []
    for _, row in df.iterrows():
        try:
            et = neptune_date_to_et(row["date"])
            sp = lro_state_moon_j2000(et)
        except spice.utils.exceptions.SpiceyError as exc:
            print(f"  SPICE error at {row['date']}: {exc}")
            continue

        sp_kep = moon_rv2coe(sp[:3], sp[3:])
        rec = dict(mjd=row["mjd"])
        for key in ("sma", "ecc", "inc", "raan", "aop", "man"):
            d = row[key] - sp_kep[key]
            if key != "sma" and key != "ecc":                     # wrap angles
                d = (d + 180.0) % 360.0 - 180.0
            rec[f"d_{key}"] = d
        records.append(rec)

    res = pd.DataFrame(records)
    print("\n--- Keplerian residuals  (NEPTUNE - SPICE/LRO) ---")
    units = dict(sma="km", ecc="-", inc="deg", raan="deg", aop="deg", man="deg")
    for key, unit in units.items():
        col = f"d_{key}"
        print(f"  {key:>4s}  RMS={_rms(res[col]):12.6f} {unit:<4s}  "
              f"max={res[col].abs().max():12.6f} {unit}")
    return res


def _rms(series):
    return float(np.sqrt((series**2).mean()))


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_mcs(res, save_path=None):
    if not HAS_MPL:
        return
    t = res["mjd"] - res["mjd"].iloc[0]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

    for col, lbl in [("dr_x", r"$\Delta r_x$"), ("dr_y", r"$\Delta r_y$"),
                     ("dr_z", r"$\Delta r_z$")]:
        ax1.plot(t, res[col] * 1e3, label=lbl)
    ax1.plot(t, res["dr_m"], color="k", lw=1.5, label=r"$|\Delta r|$")
    ax1.set_ylabel("Position residual  [m]")
    ax1.legend(ncol=4, fontsize=8); ax1.grid(True)

    for col, lbl in [("dv_x", r"$\Delta v_x$"), ("dv_y", r"$\Delta v_y$"),
                     ("dv_z", r"$\Delta v_z$")]:
        ax2.plot(t, res[col] * 1e3, label=lbl)
    ax2.plot(t, res["dv_mps"], color="k", lw=1.5, label=r"$|\Delta v|$")
    ax2.set_ylabel("Velocity residual  [mm/s]")
    ax2.set_xlabel(f"Days since MJD {res['mjd'].iloc[0]:.3f}")
    ax2.legend(ncol=4, fontsize=8); ax2.grid(True)

    fig.suptitle("NEPTUNE MCRF Cartesian vs. LRO SPICE (J2000, Moon centre)")
    fig.tight_layout()
    _save_or_show(fig, save_path)


def plot_osm(res, save_path=None):
    if not HAS_MPL:
        return
    t = res["mjd"] - res["mjd"].iloc[0]
    items = [("d_sma","SMA [km]"), ("d_ecc","ECC [-]"), ("d_inc","INC [deg]"),
             ("d_raan","RAAN [deg]"), ("d_aop","AoP [deg]"), ("d_man","Mean anom. [deg]")]

    fig, axes = plt.subplots(3, 2, figsize=(14, 9), sharex=True)
    for ax, (col, lbl) in zip(axes.flat, items):
        ax.plot(t, res[col])
        ax.axhline(0, color="k", lw=0.5)
        ax.set_ylabel(lbl); ax.grid(True)
    for ax in axes[-1]:
        ax.set_xlabel(f"Days since MJD {res['mjd'].iloc[0]:.3f}")

    fig.suptitle("NEPTUNE MCRF Keplerian elements vs. LRO SPICE")
    fig.tight_layout()
    _save_or_show(fig, save_path)


def _save_or_show(fig, path):
    if path:
        fig.savefig(path, dpi=150)
        print(f"Plot saved to {path}")
    else:
        plt.show()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--start",      default="2009-07-01T00:00:00",
                   help="LRO initial epoch  (ISO-8601 UTC, default: %(default)s)")
    p.add_argument("--duration",   type=float, default=3.0,
                   help="Propagation window in days (default: %(default)s)")
    p.add_argument("--inp",        default=str(INP_FILE),
                   help="Path to neptune.inp  (default: %(default)s)")
    p.add_argument("--mcs",        default=str(OUT_DIR / "neptune.mcs"),
                   help="NEPTUNE MCRF Cartesian output file")
    p.add_argument("--osm",        default=str(OUT_DIR / "neptune.osm"),
                   help="NEPTUNE MCRF Keplerian output file")
    p.add_argument("--no-patch",   action="store_true",
                   help="Skip patching neptune.inp  (use existing file)")
    p.add_argument("--run",        action="store_true",
                   help="Run NEPTUNE after patching neptune.inp")
    p.add_argument("--plot",       action="store_true",
                   help="Display comparison plots interactively")
    p.add_argument("--save-plots", action="store_true",
                   help="Save comparison plots as PNG files")
    args = p.parse_args()

    load_kernels()

    # ------------------------------------------------------------------
    # Step 1: extract LRO state
    # ------------------------------------------------------------------
    r_mcrf, v_mcrf, _start_et, start_dt = extract_initial_state(args.start)
    end_dt = start_dt + timedelta(days=args.duration)

    # ------------------------------------------------------------------
    # Step 2: patch neptune.inp
    # ------------------------------------------------------------------
    if not args.no_patch:
        patch_neptune_inp(args.inp, start_dt, end_dt, r_mcrf, v_mcrf)
    else:
        print("\n--no-patch set: neptune.inp not modified.")

    # ------------------------------------------------------------------
    # Step 3: run NEPTUNE
    # ------------------------------------------------------------------
    if args.run:
        run_neptune()
    else:
        print("\nNEPTUNE not run (pass --run to execute automatically).")

    # ------------------------------------------------------------------
    # Step 4: compare
    # ------------------------------------------------------------------
    ran_any = False

    if Path(args.mcs).exists():
        res_mcs = compare_mcs(args.mcs)
        if args.plot or args.save_plots:
            out = "lro_validation_mcs.png" if args.save_plots else None
            plot_mcs(res_mcs, out)
        ran_any = True
    else:
        print(f"\nMCS file not found: {args.mcs}")

    if Path(args.osm).exists():
        res_osm = compare_osm(args.osm)
        if args.plot or args.save_plots:
            out = "lro_validation_osm.png" if args.save_plots else None
            plot_osm(res_osm, out)
        ran_any = True
    else:
        print(f"OSM file not found: {args.osm}")

    if not ran_any:
        print("\nNo NEPTUNE output files found yet.  Pass --run to execute NEPTUNE,")
        print("or run it manually from the work/ directory, then re-run this script.")

    spice.kclear()


if __name__ == "__main__":
    main()
