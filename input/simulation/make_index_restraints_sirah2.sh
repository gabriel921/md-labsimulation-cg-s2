#!/bin/bash -e
# =============================================================================
# make_index_restraints_sirah2.sh — final index file + backbone restraint files
# for SIRAH 2 (the "index_final.ndx / bkbres.itp / bkbres_soft.itp" step of the
# preparation notes), including the topology edit that activates them.
#
#   1) make_ndx:  a GN GO  -> group GN_GO (backbone beads)
#   2) genrestr:  bkbres.itp (default force constants) and
#                 bkbres_soft.itp (-fc 100 100 100)
#   3) topology:  adds   #ifdef GN_GO / GN_GO_SOFT   blocks right after the
#                 "Include Position restraint file" block written by pdb2gmx.
#
#      ONE chain  : topol.top itself holds the posre.itp block
#                   -> blocks added to topol.top, files bkbres.itp / bkbres_soft.itp
#      >1 chains  : topol.top only #includes topol_<chain>.itp files, and each
#                   of those holds its own posre_<chain>.itp block
#                   -> the global restraint files are split per chain (atom
#                      numbers renumbered to be local to each chain), written
#                      as bkbres_<chain>.itp / bkbres_soft_<chain>.itp, and the
#                      blocks are added to each topol_<chain>.itp
#                      (<chain> e.g. Protein_chain_A, same as posre_<chain>.itp)
#
# Run with the current directory = 0_prep/ (the folder holding topol.top).
#
# Required environment variables:
#   GMX  GRO  TOP  NDX
# Optional:
#   BB_ATOMS   backbone bead names (default: "GN GO")
#   BB_NAME    name of the resulting index group (default: GN_GO)
# =============================================================================
: "${GMX:?GMX not set}"
: "${GRO:?GRO not set}"
: "${TOP:?TOP not set}"
: "${NDX:?NDX not set}"
BB_ATOMS="${BB_ATOMS:-GN GO}"
BB_NAME="${BB_NAME:-GN_GO}"

# -- 1) index: count the default groups so the new group's number is known ----
echo q | "$GMX" make_ndx -f "$GRO" -o index_default.ndx
N0=$(grep -c '^\[' index_default.ndx)
printf 'a %s\nname %d %s\nq\n' "$BB_ATOMS" "$N0" "$BB_NAME" | "$GMX" make_ndx -f "$GRO" -o "$NDX"
grep -q "^\[ ${BB_NAME} \]" "$NDX" || { echo "[ERROR] group ${BB_NAME} not created in $NDX" >&2; exit 1; }

# -- 2) restraint files (whole system; split per chain below if needed) -------
printf '%s\n' "$BB_NAME" | "$GMX" genrestr -f "$GRO" -n "$NDX" -o bkbres.itp
printf '%s\n' "$BB_NAME" | "$GMX" genrestr -f "$GRO" -n "$NDX" -o bkbres_soft.itp -fc 100 100 100

# -- 3) topology: include the restraint files (idempotent) --------------------

# add_restraint_blocks <topology-or-chain-itp> <hard.itp> <soft.itp>
# Inserts the GN_GO / GN_GO_SOFT blocks after the posre #endif of that file.
add_restraint_blocks() {
    local f="$1" hard="$2" soft="$3" posre_line endif_line
    if grep -q "\"$hard\"" "$f"; then
        echo "[make_index_restraints] $f already includes $hard — left unchanged."
        return 0
    fi
    posre_line=$(grep -n '#include[[:space:]]*"posre[^"]*\.itp"' "$f" | head -1 | cut -d: -f1)
    [ -n "$posre_line" ] || { echo "[ERROR] no '#include \"posre*.itp\"' block found in $f — add the GN_GO / GN_GO_SOFT blocks (see README) by hand." >&2; exit 1; }
    endif_line=$(awk -v s="$posre_line" 'NR>s && /^#endif/ {print NR; exit}' "$f")
    [ -n "$endif_line" ] || { echo "[ERROR] could not find the #endif closing the posre block in $f" >&2; exit 1; }
    {
        head -n "$endif_line" "$f"
        cat <<BLK

; Backbone restraints
#ifdef GN_GO
#include "$hard"
#endif

; Backbone soft restraints
#ifdef GN_GO_SOFT
#include "$soft"
#endif
BLK
        tail -n +"$((endif_line + 1))" "$f"
    } > "${f}.tmp" && mv "${f}.tmp" "$f"
    echo "[make_index_restraints] restraint blocks added to $f (after line $endif_line)."
}

# split_restraints <global.itp> <first_atom> <last_atom> <out.itp>
# Keeps only the restraints of atoms first..last and renumbers them from 1
# (position_restraints atom numbers are local to each [ moleculetype ]).
split_restraints() {
    awk -v s="$2" -v e="$3" '
        BEGIN { print "[ position_restraints ]"; print "; atom  type      fx      fy      fz" }
        $1 ~ /^[0-9]+$/ && $1 >= s && $1 <= e { printf "%6d %5d %8s %8s %8s\n", $1 - s + 1, $2, $3, $4, $5 }
    ' "$1" > "$4"
    grep -qE '^ *[0-9]' "$4" || { echo "[ERROR] no backbone beads found for atoms $2-$3 in $1 (-> $4)" >&2; exit 1; }
}

# Number of atoms in the [ atoms ] section of a chain .itp
count_itp_atoms() {
    awk '/^\[/ { s = ($0 ~ /^\[[[:space:]]*atoms[[:space:]]*\]/) } s && $1 ~ /^[0-9]+$/ { n++ } END { print n + 0 }' "$1"
}

TOPDIR=$(dirname "$TOP")
# chain itps in the order topol.top includes them (= order of the atoms in the .gro)
mapfile -t CHAIN_ITPS < <(grep -oE '#include[[:space:]]+"topol_[^"]+\.itp"' "$TOP" | sed -E 's/.*"([^"]+)"/\1/')

if [ "${#CHAIN_ITPS[@]}" -eq 0 ]; then
    # ---- single chain: posre block lives in topol.top itself ----------------
    add_restraint_blocks "$TOP" bkbres.itp bkbres_soft.itp
else
    # ---- several chains: one bkbres / bkbres_soft pair per chain ------------
    start=1
    for itp in "${CHAIN_ITPS[@]}"; do
        stem="${itp%.itp}"; stem="${stem#topol_}"          # e.g. Protein_chain_A
        n=$(count_itp_atoms "$TOPDIR/$itp")
        [ "$n" -gt 0 ] || { echo "[ERROR] no atoms found in [ atoms ] of $TOPDIR/$itp" >&2; exit 1; }
        end=$(( start + n - 1 ))
        echo "[make_index_restraints] $itp : atoms ${start}-${end}  -> bkbres_${stem}.itp / bkbres_soft_${stem}.itp"
        split_restraints bkbres.itp      "$start" "$end" "$TOPDIR/bkbres_${stem}.itp"
        split_restraints bkbres_soft.itp "$start" "$end" "$TOPDIR/bkbres_soft_${stem}.itp"
        add_restraint_blocks "$TOPDIR/$itp" "bkbres_${stem}.itp" "bkbres_soft_${stem}.itp"
        start=$(( end + 1 ))
    done
    rm -f bkbres.itp bkbres_soft.itp     # whole-system files are not used by any topology
fi
