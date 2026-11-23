#!/bin/bash -e
# =============================================================================
# run_setup.sh   —  SIRAH 2 coarse-grained version
# Master setup script for Aβ42 coarse-grained MD simulation (SIRAH 2)
#
# Runs sequentially:
#   0) CG mapping + topology (cgconv.pl + pdb2gmx)
#   1) Box building          (editconf)
#   2) Solvation             (solvate, SIRAH WT4 water + removal of close water)
#   3) Ionisation            (genion)  + 3b) index and backbone restraint files
#   4) Energy minimisation   (min_sirah2.sh: em1 + em2)
#   5) NVT equilibration     (eq_NVT_sirah2.sh)
#   6) NPT equilibration     (SLURM by default, or local sequential — see
#   7) Production MD          EXECUTION_MODE below)
#
# Set EXECUTION_MODE=slurm (default) or EXECUTION_MODE=local in the USER
# CONFIGURATION section to choose how steps 6-7 run. Steps 0-5 are identical
# either way.
# =============================================================================
# ── LOG EVERYTHING TO FILE ───────────────────────────────────────────────────
#LOGFILE="run_setup_$(date +%Y%m%d_%H%M%S).log"
#exec >> "$LOGFILE" 2>&1
#echo "============================================================"
#echo " Logging started: $(date)"
#echo " Log file: $LOGFILE"
#echo "============================================================"

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
PDB=ab42_1z0q.pdb                        # Input PDB file, atomistic (must be in the same directory)
SYS=ab42                                 # System name prefix for output files
N_CHAINS=1                               # Keep the first N chains from the CGCONV output PDB
TERMINAL_N=0                             # pdb2gmx -ter menu choice for the N-terminus (use the menu number)
TERMINAL_C=0                             # pdb2gmx -ter menu choice for the C-terminus (use the menu number)
GMX=/gromacs/gromacs-2018.6/bin/gmx      # gmx used for the local steps (same build as the legacy CG scripts)
export GMXLIB=$HOME/Desktop/GMXLIB       # GROMACS library path — the SIRAH force-field folder lives here

SIRAH_FF_DIR=$GMXLIB/sirah2_scaledF1_3.ff # SIRAH force-field folder (linked as sirah.ff inside every stage folder)
FF_NAME=sirah                            # pdb2gmx -ff value (the sirah.ff link created in 0_prep)
WM=1                                     # pdb2gmx water-model menu entry = "SIRAH solvent models"
PDB2GMX_EXTRA=""                         # additional pdb2gmx flags, if needed (do not put -ter here)
SOLVENT_GRO=sirah.ff/wt416.gro           # SIRAH WT4 water box (relative to 0_prep)

BOX_PADDING=3                            # Minimum distance from solute to box edge (nm)
BOX_TYPE=cubic                           # Box geometry
CLOSE_WATER_CUTOFF=0.3                   # WT4 closer than this to the protein is removed (nm)
ION_CONC=0.15                            # NaCl concentration (M)
ION_POS=NaW                              # SIRAH hydrated cation
ION_NEG=ClW                              # SIRAH hydrated anion
ION_SOLVENT_GROUP=WT4                    # group genion replaces by ions (SIRAH water)
NTHREADS=8                               # CPU threads for local steps (min, NVT)
MAXWARN=1                                # Allowed grompp warnings

MDP_DIR=$(pwd)/mdp                       # Directory containing all .mdp files
SCRIPTS_DIR=$(pwd)                       # System root — watcher script lives here
INP_DIR=$(pwd)/input/simulation          # Directory containing .dat files + slurm_*.sh + concat script

POSRES_NPT=yes                           # POSRES_NPT: backbone restraints during NPT equilibration.
                                         # yes = keep the "define" line in the NPT .mdp     no = comment it out
REPLICATES=3                             # Number of independent production replicates
                                         # Each gets its own 4_prd/repN subfolder,
                                         # all starting from the same NPT endpoint.

# ── EXECUTION MODE ────────────────────────────────────────────────────────────
# "slurm" (default) — unchanged original behaviour: NPT and production run as
#                      self-requeuing SLURM jobs (steps 6-7 below), watched in
#                      the background by watch_npt_submit_prd.sh.
# "local"            — no SLURM involved. NPT then each production replicate
#                       run sequentially, in the foreground, in this shell,
#                       reusing the same grompp/mdrun commands and the same
#                       segment/continuation logic (see run_segments_local.sh).
EXECUTION_MODE=slurm

# ── NPT / PRODUCTION STEP-COUNT AND BLOCK-COUNT CONFIGURATION ────────────────
# These control the `nsteps`/`dt` written into the copies of eq2_CGPROT.mdp
# and md_CGPROT.mdp used for the run, and the number of NPT/production
# blocks (SLURM "Segments", or local sequential iterations). Defaults below
# are a short end-to-end test (same convention as the all-atom repository) with the
# coarse-grained time step of 0.02 ps (20 fs). One segment covers nsteps x dt,
# e.g. 50000 steps x 0.02 ps = 1 ns — raise nsteps/segments for real runs.
NPT_NSTEPS=50                            # nsteps written into the NPT .mdp for this run
NPT_DT=0.02                              # dt (ps) written into the NPT .mdp for this run
NPT_SEGMENTS=1                           # number of NPT blocks (SLURM Segments= / local loop count)

PROD_NSTEPS=5000                         # nsteps written into the production .mdp for this run
PROD_DT=0.02                             # dt (ps) written into the production .mdp for this run
PROD_SEGMENTS=1                          # number of production blocks per replicate

# ── ANALYSIS INDEX GROUPS (system_groups.cfg) ────────────────────────────────
# Written once to $PREP/system_groups.cfg and sourced by every analysis
# script (index_maker.sh, Rg, Sasa, etc.) — this is the ONLY place chain and
# peptide atom ranges should need editing when you move to a different
# protein/system. See $PREP/system_groups.cfg after setup for the full
# format and comments.
# NOTE (CG): GROMACS has no default "Backbone" group for CG beads — check that your
# index_maker.sh builds the backbone group (BB for Martini 3, GN_GO for SIRAH 2).
ANALYSIS_BACKBONE_GROUP_NAME="Backbone"
# ":auto" is replaced by 1-<number of protein CG beads> once STEP 0 has run
ANALYSIS_CHAIN_DEFS=("ChainA:auto")
ANALYSIS_PEPTIDE_DEFS=()
ANALYSIS_BACKBONE_FOR=(ChainA)
ANALYSIS_ALL_CHAINS_NAME="All_Chains"

ANALYSIS_CONCAT_HAS_SOLUTE=1              # 1 = center on a solute; 0 = PBC only
ANALYSIS_CONCAT_SOLUTE_GROUP="Protein"    # group to center on (if above = 1)
ANALYSIS_CONCAT_CENTER_OUTPUT_GROUP="Protein"  # group written to the final trajectory
# ── END USER CONFIGURATION ───────────────────────────────────────────────────

# ── FOLDER STRUCTURE ─────────────────────────────────────────────────────────
PREP=0_prep
MIN=1_min
NVT=2_nvt
NPT=3_npt
PRD=4_prd

TOP=topol.top
NDX=index_final.ndx
ROOT_DIR=$(pwd)
PREP_ABS=$ROOT_DIR/$PREP

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
log()        { echo -e "\n[$(date '+%H:%M:%S')] >>> $*"; }
die()        { echo -e "\n[ERROR] $*" >&2; exit 1; }
check_file() { [ -f "$1" ] || die "Required file not found: $1"; }
check_cmd()  { command -v "$1" &>/dev/null || die "Command not found: $1 — check your GROMACS installation."; }

# Patches nsteps/dt into a copied .mdp file (never touches the template in
# $MDP_DIR — only the per-run copy passed as $1).
patch_mdp_steps() {
    local mdp_file="$1" nsteps="$2" dt="$3"
    sed -i "s/^nsteps[[:space:]]*=.*/nsteps      = ${nsteps}          ; set by run_setup.sh/" "$mdp_file"
    sed -i "s/^dt[[:space:]]*=.*/dt          = ${dt}             ; set by run_setup.sh/"     "$mdp_file"
}

# Comments out the define line in the per-run copy only (POSRES_NPT=no).
patch_mdp_posres() {
    local mdp_file="$1" use_posres="$2"
    if [ "$use_posres" = "no" ]; then
        sed -i 's/^define[[:space:]]*=.*/; define disabled by run_setup.sh (POSRES_NPT=no)/' "$mdp_file"
    fi
}

# Patches Segments= into a copied SLURM .dat file (never touches the
# template in $INP_DIR — only the per-run copy passed as $1).
patch_dat_segments() {
    local dat_file="$1" segments="$2"
    sed -i "s/^Segments=.*/Segments=${segments}/" "$dat_file"
}

# Rewrites the number of molecules <mol> in [ molecules ] (or appends the line)
# so the topology matches the .gro. Replaces the manual "grep -c <bead> file.gro
# + edit topol.top" step of the preparation notes.
set_molecule_count() {
    local top="$1" mol="$2" count="$3" tmp
    tmp="$(mktemp)"
    awk -v m="$mol" -v c="$count" '
        /^\[/ { in_mol = ($0 ~ /molecules/) }
        in_mol && $1 == m && $2 ~ /^[0-9]+$/ { printf "%-15s %5d\n", m, c; found = 1; next }
        { print }
        END { if (!found) printf "%-15s %5d\n", m, c }
    ' "$top" > "$tmp" && mv "$tmp" "$top"
}

# Number of atoms named <name> in a .gro file (atom-name field = columns 11-15).
count_atoms() {
    awk -v n="$2" 'NR > 2 { a = substr($0, 11, 5); gsub(/ /, "", a); if (a == n) c++ } END { print c + 0 }' "$1"
}

# Writes $PREP/system_groups.cfg (called right after STEP 0, once the number of
# CG beads is known — see ANALYSIS_CHAIN_DEFS above).
write_system_groups_cfg() {
    log "Writing $PREP/system_groups.cfg..."
    local def
    local -a resolved=()
    for def in "${ANALYSIS_CHAIN_DEFS[@]}"; do
        resolved+=("${def//:auto/:1-${PROTEIN_NBEADS}}")
    done
    {
        echo "# Auto-generated by run_setup.sh on $(date) — edit the ANALYSIS_*"
        echo "# arrays at the top of that script and re-run setup to regenerate."
        echo
        echo "BACKBONE_GROUP_NAME=\"${ANALYSIS_BACKBONE_GROUP_NAME}\""
        echo
        echo "CHAIN_DEFS=("
        for def in "${resolved[@]}"; do echo "    \"${def}\""; done
        echo ")"
        echo
        echo "PEPTIDE_DEFS=("
        for def in "${ANALYSIS_PEPTIDE_DEFS[@]}"; do echo "    \"${def}\""; done
        echo ")"
        echo
        echo "BACKBONE_FOR=(${ANALYSIS_BACKBONE_FOR[*]})"
        echo
        echo "ALL_CHAINS_NAME=\"${ANALYSIS_ALL_CHAINS_NAME}\""
        echo
        echo "CONCAT_HAS_SOLUTE=${ANALYSIS_CONCAT_HAS_SOLUTE}"
        echo "CONCAT_SOLUTE_GROUP=\"${ANALYSIS_CONCAT_SOLUTE_GROUP}\""
        echo "CONCAT_CENTER_OUTPUT_GROUP=\"${ANALYSIS_CONCAT_CENTER_OUTPUT_GROUP}\""
    } > "$ROOT_DIR/$PREP/system_groups.cfg"
    log "$PREP/system_groups.cfg written."
}

# Copies into a stage directory (NPT, each production replicate) everything
# grompp needs to resolve the #includes of topol.top: the restraint .itp files
# and the SIRAH links (sirah.ff, plus residuetypes.dat / specbond.dat, which
# make_ndx / trjconv need to recognise the SIRAH residues as "Protein").
stage_topology_files() {
    local dest="$1" f l
    for f in "$PREP_ABS"/*.itp; do
        if [ -e "$f" ]; then cp -f "$f" "$dest/"; fi
    done
    for l in sirah.ff residuetypes.dat specbond.dat; do
        if [ -e "$PREP_ABS/$l" ]; then ln -sfn "$(readlink -f "$PREP_ABS/$l")" "$dest/$l"; fi
    done
    return 0
}

# Checks GMXLIB (if set) and the SIRAH force-field folder. pdb2gmx itself stays
# the final authority on whether the force field is usable.
check_gmxlib() {
    if [ -n "$GMXLIB" ] && [ ! -d "$GMXLIB" ]; then
        die "GMXLIB is set to '$GMXLIB' but that directory does not exist. Fix GMXLIB in run_setup.sh."
    fi
    [ -d "$SIRAH_FF_DIR" ] || die "SIRAH force-field folder not found: $SIRAH_FF_DIR — put it inside GMXLIB or fix SIRAH_FF_DIR in run_setup.sh."
}

# ── PREFLIGHT CHECKS ─────────────────────────────────────────────────────────
log "Running preflight checks..."
check_cmd  "$GMX"
check_file "$PDB"
check_file "$MDP_DIR/em1_CGPROT.mdp"
check_file "$MDP_DIR/em2_CGPROT.mdp"
check_file "$MDP_DIR/eq1_CGPROT.mdp"
check_file "$MDP_DIR/eq2_CGPROT.mdp"
check_file "$MDP_DIR/md_CGPROT.mdp"
check_file "$INP_DIR/make_index_restraints_sirah2.sh"
check_file "$INP_DIR/min_sirah2.sh"
check_file "$INP_DIR/eq_NVT_sirah2.sh"
check_file "$INP_DIR/concat_trj_replica.sh"

case "$EXECUTION_MODE" in
    slurm)
        check_file "$INP_DIR/slurm_MD_eq_NPT_sirah2.sh"
        check_file "$INP_DIR/slurm_MD_eq_NPT_sirah2.dat"
        check_file "$INP_DIR/slurm_MD_prd_sirah2.sh"
        check_file "$INP_DIR/slurm_MD_prd_sirah2.dat"
        check_file "$SCRIPTS_DIR/watch_npt_submit_prd.sh"
        ;;
    local)
        check_file "$INP_DIR/run_segments_local.sh"
        ;;
    *)
        die "Unknown EXECUTION_MODE='$EXECUTION_MODE' — must be 'slurm' or 'local'."
        ;;
esac

case "$POSRES_NPT" in
    yes|no) ;;
    *) die "Unknown POSRES_NPT='$POSRES_NPT' — must be 'yes' or 'no'." ;;
esac

check_gmxlib
check_file "$SIRAH_FF_DIR/wt416.gro"
check_file "$SIRAH_FF_DIR/tools/CGCONV/cgconv.pl"
log "All checks passed."

# ── CREATE DIRECTORY STRUCTURE ───────────────────────────────────────────────
log "Creating directory structure..."
mkdir -p $PREP $MIN $NVT $NPT $PRD
echo "  $PREP/  — system preparation (topology, box, solvation, ions)"
echo "  $MIN/   — energy minimisation"
echo "  $NVT/   — NVT equilibration"
echo "  $NPT/   — NPT equilibration ($EXECUTION_MODE mode)"
echo "  $PRD/   — production MD ($EXECUTION_MODE mode)"

# ── STEP 0 — CG mapping (cgconv.pl) + topology (pdb2gmx) ─────────────────────
log "STEP 0 — SIRAH 2: mapping $PDB to CG and generating the topology..."
cd $PREP
ln -sfn "$SIRAH_FF_DIR"            sirah.ff            # force field, as in the notes
ln -sfn sirah.ff/residuetypes.dat  residuetypes.dat    # lets GROMACS recognise SIRAH residues as Protein
ln -sfn sirah.ff/specbond.dat      specbond.dat
./sirah.ff/tools/CGCONV/cgconv.pl -i ../$PDB -o ${SYS}_CG.pdb
check_file "${SYS}_CG.pdb"

# Keep only the first N chains in the CGCONV output (PDB chain ID column 22;
# TER records are also recognised when chain IDs are blank or reused).
case "$N_CHAINS" in
    ''|*[!0-9]*|0) die "N_CHAINS must be a positive integer." ;;
esac
awk -v n="$N_CHAINS" '
    BEGIN { chain_count=0; current_chain=""; new_chain=1 }
    /^(ATOM  |HETATM)/ {
        chain=substr($0,22,1)
        if (new_chain || chain != current_chain) {
            chain_count++
            current_chain=chain
            new_chain=0
        }
        if (chain_count > n) exit
        print
        next
    }
    /^TER/ {
        if (chain_count <= n) print
        new_chain=1
        current_chain=""
        next
    }
    chain_count <= n { print }
' "${SYS}_CG.pdb" > "${SYS}_CG_selected.pdb"
mv "${SYS}_CG_selected.pdb" "${SYS}_CG.pdb"
CG_CHAIN_COUNT=$(awk '
    BEGIN { count=0; current_chain=""; new_chain=1 }
    /^(ATOM  |HETATM)/ {
        chain=substr($0,22,1)
        if (new_chain || chain != current_chain) {
            count++
            current_chain=chain
            new_chain=0
        }
    }
    /^TER/ { new_chain=1; current_chain="" }
    END { print count+0 }
' "${SYS}_CG.pdb")
[ "$CG_CHAIN_COUNT" -eq "$N_CHAINS" ] || die "Requested $N_CHAINS chain(s), but only $CG_CHAIN_COUNT chain(s) were found in ${SYS}_CG.pdb. Check the PDB chain IDs/TER records."
log "Keeping the first $CG_CHAIN_COUNT chain(s) from ${SYS}_CG.pdb."

# pdb2gmx asks for the water model first, then the N- and C-terminal
# choices for each chain when -ter is enabled.
case "$TERMINAL_N:$TERMINAL_C" in
    *[!0-9:]*|:*) die "TERMINAL_N and TERMINAL_C must be non-negative menu numbers." ;;
esac
PDB2GMX_INPUT=("$WM")
for ((chain=1; chain<=CG_CHAIN_COUNT; chain++)); do
    PDB2GMX_INPUT+=("$TERMINAL_N" "$TERMINAL_C")
done
printf '%s\n' "${PDB2GMX_INPUT[@]}" | \
    $GMX pdb2gmx -f ${SYS}_CG.pdb -o ${SYS}_CG.gro -ff $FF_NAME -ter $PDB2GMX_EXTRA
check_file "${SYS}_CG.gro"
check_file "$TOP"
grep -qi sirah "$TOP" || die "topol.top does not reference the SIRAH force field — check FF_NAME/WM and the sirah.ff link in $PREP/."
PROTEIN_NBEADS=$(sed -n 2p "${SYS}_CG.gro" | tr -d '[:space:]')
log "Protein CG beads: $PROTEIN_NBEADS"
write_system_groups_cfg

# ── STEP 1 — Build simulation box ────────────────────────────────────────────
log "STEP 1 — editconf: building $BOX_TYPE box (padding = ${BOX_PADDING} nm)..."
$GMX editconf \
    -f  ${SYS}_CG.gro \
    -o  ${SYS}_box.gro \
    -c \
    -d  $BOX_PADDING \
    -bt $BOX_TYPE
check_file "${SYS}_box.gro"

# ── STEP 2 — Solvation (WT4) + removal of water too close to the protein ─────
log "STEP 2 — solvate: filling box with SIRAH WT4 water..."
$GMX solvate \
    -cp ${SYS}_box.gro \
    -cs $SOLVENT_GRO \
    -o  ${SYS}_solv.gro
check_file "${SYS}_solv.gro"
set_molecule_count "$TOP" WT4 "$(count_atoms "${SYS}_solv.gro" WP1)"

log "STEP 2b — removing WT4 within ${CLOSE_WATER_CUTOFF} nm of the protein..."
echo q | $GMX make_ndx -f ${SYS}_solv.gro -o index_tmp.ndx
$GMX grompp \
    -f $MDP_DIR/em1_CGPROT.mdp \
    -p $TOP \
    -po delete1.mdp \
    -c ${SYS}_solv.gro \
    -o solv.tpr \
    -maxwarn $MAXWARN
$GMX select \
    -f ${SYS}_solv.gro \
    -s solv.tpr \
    -n index_tmp.ndx \
    -on rm_close_wt4.ndx \
    -select "not (same residue as (resname WT4 and within ${CLOSE_WATER_CUTOFF} of group Protein))"
echo 0 | $GMX editconf \
    -f ${SYS}_solv.gro \
    -o ${SYS}_solv_clean.gro \
    -n rm_close_wt4.ndx
check_file "${SYS}_solv_clean.gro"
set_molecule_count "$TOP" WT4 "$(count_atoms "${SYS}_solv_clean.gro" WP1)"

# ── STEP 3 — Add ions ────────────────────────────────────────────────────────
log "STEP 3 — genion: neutralising system and adding ${ION_CONC} M NaCl..."
$GMX grompp \
    -f $MDP_DIR/em1_CGPROT.mdp \
    -p $TOP \
    -po delete2.mdp \
    -c ${SYS}_solv_clean.gro \
    -o solvfinal.tpr \
    -maxwarn $MAXWARN

echo "$ION_SOLVENT_GROUP" | $GMX genion \
    -s solvfinal.tpr \
    -o ${SYS}_solv_ions.gro \
    -p $TOP \
    -pname $ION_POS \
    -nname $ION_NEG \
    -neutral \
    -conc $ION_CONC
check_file "${SYS}_solv_ions.gro"
# make [ molecules ] match the .gro exactly (the notes edit these counts by hand)
set_molecule_count "$TOP" WT4      "$(count_atoms "${SYS}_solv_ions.gro" WP1)"
set_molecule_count "$TOP" $ION_POS "$(count_atoms "${SYS}_solv_ions.gro" $ION_POS)"
set_molecule_count "$TOP" $ION_NEG "$(count_atoms "${SYS}_solv_ions.gro" $ION_NEG)"

# ── STEP 3b — Index file + backbone restraint files ──────────────────────────
log "STEP 3b — index file and backbone restraint files (make_ndx / genrestr)..."
GMX=$GMX \
GRO=${SYS}_solv_ions.gro \
TOP=$TOP \
NDX=$NDX \
    bash $INP_DIR/make_index_restraints_sirah2.sh
check_file "$NDX"
compgen -G "bkbres*.itp"      >/dev/null || die "No bkbres*.itp generated in $PREP/"
compgen -G "bkbres_soft*.itp" >/dev/null || die "No bkbres_soft*.itp generated in $PREP/"
cd ..

# ── STEP 4 — Energy minimisation ─────────────────────────────────────────────
log "STEP 4 — Energy minimisation..."
cd $MIN
GMX=$GMX \
NTHREADS=$NTHREADS \
MDP=$MDP_DIR/em1_CGPROT.mdp \
MDP2=$MDP_DIR/em2_CGPROT.mdp \
GRO=../$PREP/${SYS}_solv_ions.gro \
TOP=../$PREP/$TOP \
NDX=../$PREP/$NDX \
    bash $INP_DIR/min_sirah2.sh
check_file "em.gro"
cd ..

# ── STEP 5 — NVT equilibration ───────────────────────────────────────────────
log "STEP 5 — NVT equilibration..."
cd $NVT
cp -f ../$PREP/$NDX .          # index created in STEP 3b
GMX=$GMX \
NTHREADS=$NTHREADS \
MDP=$MDP_DIR/eq1_CGPROT.mdp \
GRO=../$MIN/em.gro \
TOP=../$PREP/$TOP \
NDX=$NDX \
    bash $INP_DIR/eq_NVT_sirah2.sh
check_file "nvt.gro"
cd ..

if [ "$EXECUTION_MODE" = "slurm" ]; then
    # ── STEP 6 — NPT equilibration (SLURM) ───────────────────────────────────
    log "STEP 6 — Setting up NPT folder and submitting SLURM job..."
    cd $NPT
    cp ../$PREP/$TOP       .
    stage_topology_files .
    cp ../$NVT/$NDX        .
    cp ../$NVT/nvt.gro     .
    cp ../$NVT/nvt.cpt     .
    cp $MDP_DIR/eq2_CGPROT.mdp .
    patch_mdp_steps eq2_CGPROT.mdp "$NPT_NSTEPS" "$NPT_DT"
    patch_mdp_posres eq2_CGPROT.mdp "$POSRES_NPT"
    cp $INP_DIR/slurm_MD_eq_NPT_sirah2.dat .
    patch_dat_segments slurm_MD_eq_NPT_sirah2.dat "$NPT_SEGMENTS"
    cp $INP_DIR/slurm_MD_eq_NPT_sirah2.sh  .
    chmod +x slurm_MD_eq_NPT_sirah2.sh
    ./slurm_MD_eq_NPT_sirah2.sh
    log "NPT SLURM job submitted from $(pwd)."
    cd ..

    # ── STEP 7 — Production MD folders (one per replicate) + background watcher ──
    log "STEP 7 — Preparing $REPLICATES production replicate folder(s)..."
    cd $PRD

    PRD_REP_DIRS=()   # absolute paths, filled in below, passed to the watcher

    for r in $(seq -w 1 "$REPLICATES"); do
        REPDIR="rep${r}"
        mkdir -p "$REPDIR"
        cp ../$PREP/$TOP               "$REPDIR/"
        stage_topology_files "$REPDIR"
        cp ../$NPT/$NDX                "$REPDIR/"
        cp $MDP_DIR/md_CGPROT.mdp           "$REPDIR/"
        patch_mdp_steps "$REPDIR/md_CGPROT.mdp" "$PROD_NSTEPS" "$PROD_DT"
        cp $INP_DIR/slurm_MD_prd_sirah2.dat  "$REPDIR/"
        patch_dat_segments "$REPDIR/slurm_MD_prd_sirah2.dat" "$PROD_SEGMENTS"
        cp $INP_DIR/slurm_MD_prd_sirah2.sh   "$REPDIR/"
        cp $INP_DIR/concat_trj_replica.sh     "$REPDIR/"

        # Suffix the job Name so each replicate shows up distinctly in squeue,
        # e.g. Name=ab42_prd  ->  Name=ab42_prd_rep1
        if grep -qE '^[[:space:]]*Name=' "$REPDIR/slurm_MD_prd_sirah2.dat"; then
            sed -i -E "s|^([[:space:]]*Name=[^[:space:]#]+)|\1_rep${r}|" "$REPDIR/slurm_MD_prd_sirah2.dat"
        else
            log "WARNING: no Name= line found in $REPDIR/slurm_MD_prd_sirah2.dat — job names will collide in squeue."
        fi

        chmod +x "$REPDIR/slurm_MD_prd_sirah2.sh"
        PRD_REP_DIRS+=("$(pwd)/$REPDIR")
        echo "  $PRD/$REPDIR/  — replicate $r"
    done
    log "$REPLICATES production folder(s) ready."
    cd ..

    # Launch a single watcher in the background, pointed at the absolute NPT path
    # and every replicate PRD path, so it survives any later cd changes in the
    # calling shell. It waits once for NPT to finish, then fans out production
    # submissions across all replicates, checking each off as it submits.
    # NPT_BASE / PRD_BASE tell the (shared) watcher which file names to look
    # for — the all-atom defaults are slurm_MD_*_allatom.*, so SIRAH overrides.
    WATCHER_LOG="watch_npt_submit_prd_$(date +%Y%m%d_%H%M%S).log"
    NPT_BASE=slurm_MD_eq_NPT_sirah2 PRD_BASE=slurm_MD_prd_sirah2 \
    nohup bash "$SCRIPTS_DIR/watch_npt_submit_prd.sh" \
          "$(pwd)/$NPT" \
          "${PRD_REP_DIRS[@]}" \
          >> "$WATCHER_LOG" 2>&1 &
    WATCHER_PID=$!

    # ── SUMMARY ──────────────────────────────────────────────────────────────────
    echo ""
    echo "============================================================"
    echo "  Setup complete — $(date)"
    echo "============================================================"
    echo ""
    echo "  Directory layout:"
    echo "    $PREP/  — topology, solvated+ionised structure"
    echo "    $MIN/   — energy minimisation output (em.gro)"
    echo "    $NVT/   — NVT equilibration output  (nvt.gro)"
    echo "    $NPT/   — NPT equilibration         (SLURM job submitted, $NPT_SEGMENTS segment(s))"
    echo "    $PRD/   — production MD             ($REPLICATES replicate(s), auto-submitted, $PROD_SEGMENTS segment(s) each)"
    for r in $(seq -w 1 "$REPLICATES"); do
        echo "      $PRD/rep${r}/"
    done
    echo ""
    echo "  ✔  Background watcher started (PID $WATCHER_PID)"
    echo "     It will poll $NPT/ every 20 s for the last NPT segment .xtc,"
    echo "     then submit production for all $REPLICATES replicate(s):"
    for r in $(seq -w 1 "$REPLICATES"); do
        echo "       - $PRD/rep${r}/slurm_MD_prd_sirah2.sh"
    done
    echo "     Once each replicate's production finishes, the same watcher"
    echo "     automatically concatenates its trajectory (concat_trj_replica.sh)"
    echo "     into $PRD/repN/traj_all_center_PBC.xtc — the file the analysis"
    echo "     scripts (Rg, etc.) expect to find."
    echo "     Watcher log: $WATCHER_LOG"
    echo ""
    echo "  To monitor the watcher:"
    echo "     tail -f $WATCHER_LOG"
    echo ""
    echo "  To cancel the watcher (if needed):"
    echo "     kill $WATCHER_PID"
    echo ""

else
    # ── LOCAL MODE — STEP 6 — NPT equilibration (sequential, foreground) ────────
    log "STEP 6 — NPT equilibration (LOCAL mode, $NPT_SEGMENTS segment(s))..."
    cd $NPT
    cp ../$PREP/$TOP       .
    stage_topology_files .
    cp ../$NVT/$NDX        .
    cp ../$NVT/nvt.gro     .
    cp ../$NVT/nvt.cpt     .
    cp $MDP_DIR/eq2_CGPROT.mdp .
    patch_mdp_steps eq2_CGPROT.mdp "$NPT_NSTEPS" "$NPT_DT"
    patch_mdp_posres eq2_CGPROT.mdp "$POSRES_NPT"

    GMX=$GMX NTHREADS=$NTHREADS MDP=eq2_CGPROT.mdp TOP=$TOP NDX=$NDX \
    GRO=nvt.gro SEGMENTS=$NPT_SEGMENTS MAXWARN=$MAXWARN REF_PREV=1 \
        bash $INP_DIR/run_segments_local.sh

    NPT_LAST_SEG=$(printf "%03d" "$NPT_SEGMENTS")
    check_file "${NPT_LAST_SEG}.gro"
    NPT_FINAL_GRO="$(pwd)/${NPT_LAST_SEG}.gro"
    log "NPT complete (local): $NPT_FINAL_GRO"
    cd ..

    # ── LOCAL MODE — STEP 7 — Production MD, replicates run one after another ──
    log "STEP 7 — Production MD (LOCAL mode, $REPLICATES replicate(s), $PROD_SEGMENTS segment(s) each)..."
    cd $PRD

    for r in $(seq -w 1 "$REPLICATES"); do
        REPDIR="rep${r}"
        mkdir -p "$REPDIR"
        cp ../$PREP/$TOP        "$REPDIR/"
        stage_topology_files "$REPDIR"
        cp ../$NPT/$NDX         "$REPDIR/"
        cp $MDP_DIR/md_CGPROT.mdp "$REPDIR/"
        patch_mdp_steps "$REPDIR/md_CGPROT.mdp" "$PROD_NSTEPS" "$PROD_DT"
        cp $INP_DIR/concat_trj_replica.sh "$REPDIR/"
        cp "$NPT_FINAL_GRO" "$REPDIR/npt_final.gro"

        log "  Replicate $r — running $PROD_SEGMENTS segment(s) sequentially..."
        (
            cd "$REPDIR"
            GMX=$GMX NTHREADS=$NTHREADS MDP=md_CGPROT.mdp TOP=$TOP NDX=$NDX \
            GRO=npt_final.gro SEGMENTS=$PROD_SEGMENTS MAXWARN=$MAXWARN \
                bash $INP_DIR/run_segments_local.sh
        )

        log "  Replicate $r — concatenating trajectory..."
        ( cd "$REPDIR" && GROM=$GMX bash concat_trj_replica.sh )

        echo "  $PRD/$REPDIR/  — replicate $r complete (local)"
    done
    cd ..

    # ── SUMMARY ──────────────────────────────────────────────────────────────
    echo ""
    echo "============================================================"
    echo "  Setup complete — $(date)  [LOCAL execution mode]"
    echo "============================================================"
    echo ""
    echo "  Directory layout:"
    echo "    $PREP/  — topology, solvated+ionised structure"
    echo "    $MIN/   — energy minimisation output (em.gro)"
    echo "    $NVT/   — NVT equilibration output  (nvt.gro)"
    echo "    $NPT/   — NPT equilibration         ($NPT_SEGMENTS segment(s), completed)"
    echo "    $PRD/   — production MD             ($REPLICATES replicate(s), $PROD_SEGMENTS segment(s) each, completed)"
    for r in $(seq -w 1 "$REPLICATES"); do
        echo "      $PRD/rep${r}/traj_all_center_PBC.xtc  — ready for run_setup_analysis.sh"
    done
    echo ""
    echo "  All stages ran sequentially in this shell — no SLURM was used."
    echo ""
fi
