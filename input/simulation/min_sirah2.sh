#!/bin/bash -e
# =============================================================================
# min_sirah2.sh — two-stage energy minimisation of the SIRAH 2 system, exactly
# the sequence in the preparation notes:
#   em1  minimisation with the backbone restrained   (em1_CGPROT.mdp, -r, -n)
#   em2  minimisation of the entire system           (em2_CGPROT.mdp)
#
# Driven by environment variables like the all-atom min script so run_setup.sh
# can call it. Run with the current directory = 1_min/.
#
# Required environment variables:
#   GMX  NTHREADS  MDP (=em1)  MDP2 (=em2)  GRO  TOP  NDX
# Output: em.gro (final structure) — min_CG.gro is the intermediate em1 result.
# =============================================================================
: "${GMX:?GMX not set}"
: "${NTHREADS:?NTHREADS not set}"
: "${MDP:?MDP (em1_CGPROT.mdp) not set}"
: "${MDP2:?MDP2 (em2_CGPROT.mdp) not set}"
: "${GRO:?GRO not set}"
: "${TOP:?TOP not set}"
: "${NDX:?NDX not set}"

# em1 — side chains/solvent relax with the backbone restrained
"$GMX" grompp -f "$MDP" -p "$TOP" -po em1.mdp -n "$NDX" -c "$GRO" -r "$GRO" -o min_CG.tpr
"$GMX" mdrun -nt "$NTHREADS" -deffnm min_CG

# em2 — whole system
"$GMX" grompp -f "$MDP2" -p "$TOP" -po em2.mdp -n "$NDX" -c min_CG.gro -o em.tpr
"$GMX" mdrun -nt "$NTHREADS" -deffnm em
