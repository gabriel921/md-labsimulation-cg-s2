#!/bin/bash -e
# =============================================================================
# script_SASA_replica.sh
#
# SIRAH 2 version. Runs gmx sasa on the "Protein" group for the single
# replica this script is run from, then produces a SASA vs time plot +
# probability-density distribution via plot_sasa.py.
#
# SIRAH difference: GROMACS has no van der Waals radii for SIRAH bead names,
# so a vdwradii.dat is generated here (one line per bead name in the .gro,
# all with the same radius, BEAD_RADIUS) and gmx reads it from the current
# directory. This is a COARSE, UNIFORM-RADIUS estimate: use it to compare
# replicas/conditions, not as an absolute SASA. For absolute values, compute
# SASA on backmapped structures with the all-atom script (see README).
#
# Run from inside a replica's analysis subfolder, e.g. 4_prd/rep1/analysis/05_SASA/
#
# Env vars (all optional, fallback defaults match that layout):
#   GROM       - gmx binary                 (default: gmx)
#   PRD_DIR    - replica root directory     (default: ../../)
#   TPR        - .tpr file                  (default: auto-detected in PRD_DIR)
#   TRAJ       - trajectory file            (default: auto-detected in PRD_DIR)
#   REP_LABEL  - label for output filenames (default: basename of PRD_DIR)
#   OUT_DIR    - where outputs go           (default: current directory)
#   BIN_WIDTH  - override auto (Freedman-Diaconis) distribution bin width
#   GRO        - structure used to list bead names (default: first NNN.gro in PRD_DIR)
#   BEAD_RADIUS- radius (nm) given to every bead (default: 0.21)
#   PROBE      - probe radius (nm)          (default: 0.21 = WT4 bead radius)
# =============================================================================

set -e

GROM="${GROM:-gmx}"
NDX="${NDX:-../00_index/index_analysis.ndx}"
PRD_DIR="${PRD_DIR:-../..}"
OUT_DIR="${OUT_DIR:-.}"
GROUP="Protein"
BEAD_RADIUS="${BEAD_RADIUS:-0.21}"
PROBE="${PROBE:-0.21}"

log() { echo -e "[$(date '+%H:%M:%S')] $*"; }
die() { echo -e "[ERROR] $*" >&2; exit 1; }

mkdir -p "$OUT_DIR"
[ -d "$PRD_DIR" ] || die "replica directory not found: $PRD_DIR"
REP_LABEL="${REP_LABEL:-$(basename "$(cd "$PRD_DIR" && pwd)")}"

# ── AUTO-DETECT TPR ────────────────────────────────────────────────────────────
if [ -z "$TPR" ]; then
    TPR=$(find "$PRD_DIR" -maxdepth 1 -name "*.tpr" | head -n1)
fi
[ -n "$TPR" ] && [ -f "$TPR" ] || die "no .tpr found in $PRD_DIR — set TPR explicitly."

# ── LOCATE TRAJECTORY ──────────────────────────────────────────────────────────
# Same fixed name concat_trj_replica.sh always writes — no fallback to a raw
# segment .xtc, since that would silently run SASA on incomplete data with
# no warning (same convention as script_Rg_replica.sh).
if [ -z "$TRAJ" ]; then
    TRAJ="$PRD_DIR/traj_all_center_PBC.xtc"
fi
[ -f "$TRAJ" ] || die "concatenated trajectory not found: $TRAJ — has concat_trj_replica.sh run for this replica yet? (check for a .concat_done marker in $PRD_DIR, or set TRAJ explicitly)"

[ -f "$NDX" ] || die "index file not found: $NDX (run index_maker.sh first)"

if [ -z "$GRO" ]; then
    GRO=$(find "$PRD_DIR" -maxdepth 1 -name '[0-9]*.gro' | sort | head -n1)
fi
[ -n "$GRO" ] && [ -f "$GRO" ] || die "no segment .gro found in $PRD_DIR — set GRO explicitly."

# ── SIRAH BEAD RADII (read by gmx sasa from the current directory) ────────────
awk -v r="$BEAD_RADIUS" 'NR > 2 && NF >= 4 { a = substr($0, 11, 5); gsub(/ /, "", a); if (a != "" && !(a in s)) { s[a] = 1; printf "??? %s %s\n", a, r } }' "$GRO" > vdwradii.dat
[ -s vdwradii.dat ] || die "could not build vdwradii.dat from $GRO"

log "==============================================="
log "SASA analysis — replica: $REP_LABEL"
log "==============================================="
log "TPR   : $TPR"
log "Traj  : $TRAJ"
log "Index : $NDX"
log "Group : $GROUP"
log "Radii : $BEAD_RADIUS nm per bead ($(wc -l < vdwradii.dat) bead names, ./vdwradii.dat), probe $PROBE nm"
echo

outfile="$OUT_DIR/SASA_${GROUP}_${REP_LABEL}.xvg"
log "→ group '$GROUP' → $outfile"
echo "$GROUP" | "$GROM" sasa \
    -s "$TPR" \
    -n "$NDX" \
    -f "$TRAJ" \
    -probe "$PROBE" \
    -o "$outfile"

log "✅ Finished SASA calculation for replica $REP_LABEL"

# ── PLOT + DISTRIBUTION ────────────────────────────────────────────────────────
PLOT_SCRIPT="$(dirname "$0")/plot_sasa.py"
if [ -f "$PLOT_SCRIPT" ]; then
    echo
    log "Generating SASA vs time + distribution plots..."
    BIN_ARGS=()
    [ -n "$BIN_WIDTH" ] && BIN_ARGS=(--bin-width "$BIN_WIDTH")
    python3 "$PLOT_SCRIPT" --outdir "$OUT_DIR" --label "$REP_LABEL" "${BIN_ARGS[@]}" "$outfile"
else
    log "WARNING: plot_sasa.py not found next to this script — skipping plot."
fi