# Aβ42 Coarse-Grained MD Workflow (SIRAH 2, GROMACS)

End-to-end coarse-grained molecular dynamics workflow using the SIRAH 2 force
field: starts from an atomistic PDB structure, maps it to CG, runs system
preparation, energy minimisation, NVT/NPT equilibration and production MD,
then hands off to a separate analysis pipeline. This script works for Aβ42;
to change the system just use a PDB structure from a different protein.
`run_setup.sh` counts the CG beads automatically after the mapping, so the
chain definition (`ANALYSIS_CHAIN_DEFS=("ChainA:auto")`) only needs editing
for multi-chain systems (set `N_CHAINS` accordingly — see
[Chains and termini](#chains-and-termini-n_chains-termini_n-termini_c)).

Two scripts drive everything:

- **`run_setup.sh`** — master setup + simulation script (CG mapping → prep →
  min → NVT → NPT → production).
- **`run_setup_analysis.sh`** — master analysis script, run once production
  trajectories are ready.

The simulation stages (NPT + production) can run either as **SLURM jobs**
(the original, cluster-based workflow) or **locally**, sequentially, in your
own shell, with no SLURM involved. Both modes reuse the exact same
`grompp`/`mdrun` commands and the same segment-continuation logic — only the
submission/dependency mechanism differs.

Companion repository: the all-atom (CHARMM36m) version of this workflow,
with the same structure and configuration logic.

---

## 1. Requirements

- GROMACS (tested with 2018.6; both CPU and GPU builds are supported)
- Perl (used by the SIRAH `cgconv.pl` mapping tool)
- Bash
- The SIRAH 2 force-field folder (see [Force fields and GMXLIB](#6-force-fields-and-gmxlib))
- For SLURM mode: a SLURM cluster with `sbatch`, and either a CPU or GPU
  partition set up (see [SLURM execution](#5-slurm-execution))
- Optional: VMD (+ AmberTools for backmapping) — see
  [SIRAH Tools in VMD](#8-sirah-tools-in-vmd)
- The analysis scripts have their own dependencies — see the comments at the
  top of `run_setup_analysis.sh` and the scripts it calls.

---

## 2. Installation

Clone the repository and make sure the expected layout is present:

```
.
├── run_setup.sh
├── run_setup_analysis.sh
├── watch_npt_submit_prd.sh
├── mdp/
│   ├── em1_CGPROT.mdp
│   ├── em2_CGPROT.mdp
│   ├── eq1_CGPROT.mdp
│   ├── eq2_CGPROT.mdp
│   └── md_CGPROT.mdp
└── input/
    ├── simulation/
    │   ├── make_index_restraints_sirah2.sh
    │   ├── min_sirah2.sh
    │   ├── eq_NVT_sirah2.sh
    │   ├── run_segments_local.sh
    │   ├── concat_trj_replica.sh
    │   ├── slurm_MD_eq_NPT_sirah2.sh
    │   ├── slurm_MD_eq_NPT_sirah2.dat
    │   ├── slurm_MD_prd_sirah2.sh
    │   └── slurm_MD_prd_sirah2.dat
    └── analysis/
        ├── 00_index/index_maker.sh
        ├── 01_Rg/ ...
        └── (one subfolder per analysis — see run_setup_analysis.sh)
```

Place your input (atomistic) PDB file in the repository root (default filename
`ab42_1z0q.pdb`, configurable — see below).

No build step is required; the scripts call your GROMACS installation
directly.

---

## 3. Configuration

All user-facing configuration lives in the **USER CONFIGURATION** block at
the top of `run_setup.sh`. Nothing below that block needs editing for
normal use.

### Core settings

| Variable | Meaning |
|---|---|
| `PDB` | Input atomistic PDB filename (must sit next to the script) |
| `SYS` | System name prefix used for output files |
| `GMX` | Path to the `gmx` executable used for local steps (prep, minimisation, NVT, and — in local mode — NPT/production) |
| `GMXLIB` | Force-field library directory — see [Force fields and GMXLIB](#6-force-fields-and-gmxlib) |
| `SIRAH_FF_DIR` | SIRAH 2 force-field folder (linked as `sirah.ff` inside every stage folder) |
| `FF_NAME`, `WM` | Force-field name and water-model menu entry passed to `pdb2gmx` |
| `PDB2GMX_EXTRA` | Extra `pdb2gmx` flags (do **not** put `-ter` here, it is already set) |
| `SOLVENT_GRO` | SIRAH WT4 water box used for solvation |
| `BOX_PADDING`, `BOX_TYPE` | Box parameters |
| `CLOSE_WATER_CUTOFF` | WT4 closer than this to the protein is removed (nm) |
| `ION_CONC`, `ION_POS`, `ION_NEG`, `ION_SOLVENT_GROUP` | Ionisation parameters (SIRAH hydrated ions, replacing WT4) |
| `NTHREADS` | Thread count for local steps |
| `MAXWARN` | `grompp -maxwarn` value |
| `REPLICATES` | Number of independent production replicates (each gets its own `4_prd/repN/`) |

### Chains and termini (`N_CHAINS`, `TERMINAL_N`, `TERMINAL_C`)

`pdb2gmx` is run with the `-ter` flag, so it asks, for every chain, which
start (N-terminus) and end (C-terminus) type to use. The script answers these
prompts for you, based on three variables in the USER CONFIGURATION block:

```bash
N_CHAINS=1       # chains kept from the CGCONV output PDB
TERMINAL_N=0     # start (N-terminus) menu number
TERMINAL_C=0     # end (C-terminus) menu number
```

`N_CHAINS` is the number of chains to keep from the mapped CG PDB (the first
N, by chain ID / `TER` records); the script stops if fewer are found. The
same `TERMINAL_N` / `TERMINAL_C` answers are repeated for every chain, so all
chains get the same terminus choice.

> **Warning.** The menu numbers depend on the SIRAH force-field version.
> Always check manually before trusting these settings: run
> `gmx pdb2gmx -f <your_CG.pdb> -ff sirah -ter` by hand, read the menu for your
> own system, and confirm which numbers correspond to the charged/neutral
> termini.

### Position restraints in NPT (`POSRES_NPT`)

```bash
POSRES_NPT=yes   # or: no
```

- `yes` — keeps the `define` line in the NPT `.mdp`, so the protein backbone
  is restrained (`bkbres.itp`) during NPT equilibration.
- `no` — comments that line out, so NPT runs without restraints.

Only the per-run copy of `eq2_CGPROT.mdp` placed in `3_npt/` is modified;
the template in `mdp/` is never touched. With `POSRES_NPT=yes`, the `grompp`
call in the NPT scripts needs a reference structure for the restraints
(`-r`), e.g. the NVT output `nvt.gro`. For multi-chain systems, one
`bkbres` / `bkbres_soft` file pair is generated per chain.

### Execution mode

```bash
EXECUTION_MODE=slurm   # or: local
```

- `slurm` (default) — NPT and production run as self-requeuing SLURM jobs,
  with a background watcher chaining them.
- `local` — no SLURM: NPT then each production replicate run sequentially,
  in the foreground, in the current shell.

See [Local execution](#4-local-execution) and
[SLURM execution](#5-slurm-execution) below.

### Step counts and block counts

```bash
NPT_NSTEPS=50           # nsteps written into the NPT .mdp for this run
NPT_DT=0.02             # dt (ps) written into the NPT .mdp for this run
NPT_SEGMENTS=1          # number of NPT blocks (SLURM Segments= / local loop count)

PROD_NSTEPS=5000        # nsteps written into the production .mdp for this run
PROD_DT=0.02            # dt (ps) written into the production .mdp for this run
PROD_SEGMENTS=1         # number of production blocks per replicate
```

These control:

- The `nsteps` and `dt` values written into the **copies** of
  `eq2_CGPROT.mdp` / `md_CGPROT.mdp` used for the run (the template files
  under `mdp/` are never modified).
- The number of NPT/production blocks — the SLURM `Segments=` value in
  `slurm_MD_eq_NPT_sirah2.dat` / `slurm_MD_prd_sirah2.dat`, or the number
  of sequential loop iterations in local mode.

One block/segment covers `nsteps × dt` of simulation time and continues
directly from the previous block's output structure — this is unchanged
in both execution modes. The defaults above are a short end-to-end test with
the CG time step of 0.02 ps (20 fs); for real runs raise `nsteps`/segments
(e.g. 50000 steps × 0.02 ps = 1 ns per segment).

---

## 4. Local execution

Set `EXECUTION_MODE=local` and run:

```bash
bash run_setup.sh
```

Steps 0–5 (CG mapping and topology, box, solvation, ionisation, minimisation,
NVT) are identical to SLURM mode. For steps 6–7, the script calls
`input/simulation/run_segments_local.sh` directly — once for NPT, then once
per production replicate, one after another, entirely in the foreground.
Nothing is submitted to a queue and no background watcher is started.

`run_segments_local.sh` reimplements the same per-segment logic the SLURM
scripts use (segment numbering from existing `NNN.gro` files, `tinit`
computed as `nsteps × dt × (segment − 1)` read straight from the `.mdp`,
`grompp` + `mdrun` continuing from the previous segment) — just without
`sbatch`, scratch staging, or self-requeue. It is safe to re-run: any
segment whose `.gro` already exists is skipped.

Once production finishes, each replicate's trajectory is concatenated
automatically (`concat_trj_replica.sh`), producing
`4_prd/repN/traj_all_center_PBC.xtc` — the same file the analysis scripts
expect in SLURM mode.

Local mode is intended for workstations or interactive nodes; it uses the
`GMX`/`NTHREADS` values from the configuration block for every stage,
including NPT and production.

---

## 5. SLURM execution

Set `EXECUTION_MODE=slurm` (the default) and run:

```bash
./run_setup.sh
```

1. NPT equilibration is submitted via `slurm_MD_eq_NPT_sirah2.sh`, which
   self-requeues until `Segments=` (from `slurm_MD_eq_NPT_sirah2.dat`) is
   reached.
2. `watch_npt_submit_prd.sh` runs in the background, waits for NPT's last
   segment, then submits production for every replicate
   (`slurm_MD_prd_sirah2.sh`, self-requeuing the same way).
3. Once each replicate's production finishes, the watcher concatenates its
   trajectory automatically.

The step-count/block-count configuration variables described above are
applied to the copies of the `.mdp`/`.dat` files placed in `3_npt/` and each
`4_prd/repN/` *before* the SLURM scripts are invoked, so the dependency chain,
self-requeue logic, and directory layout are completely unaffected.

### Running NPT equilibration on GPU

To run the NPT equilibration on a GPU, edit
`input/simulation/slurm_MD_eq_NPT_sirah2.dat` and change both the GROMACS
version (to a GPU-enabled build) and the partition (to a GPU partition on your
cluster). Nothing else in the workflow needs to change.

Monitor with:

```bash
tail -f watch_npt_submit_prd_*.log
squeue -u $USER
```

---

## 6. Force fields and GMXLIB

SIRAH is not bundled with GROMACS. `GMXLIB` points GROMACS at an external
directory that holds the SIRAH 2 force-field folder:

```bash
export GMXLIB="/path/to/my/gromacs-topologies"
SIRAH_FF_DIR=$GMXLIB/sirah2.ff     # your SIRAH 2 folder name
```

`run_setup.sh` links `SIRAH_FF_DIR` as `sirah.ff` inside every stage folder
(`0_prep/`, `3_npt/`, each `4_prd/repN/`), together with `residuetypes.dat` and
`specbond.dat`, so GROMACS tools recognise the SIRAH residues as Protein.
`run_setup_analysis.sh` makes the same links in every analysis folder.

**Notes:**

- `GMXLIB` should point to the directory *containing* the force-field folder.
- The force-field folder must keep the internal structure SIRAH ships with
  (`forcefield.itp`, `wt416.gro`, `tools/CGCONV/cgconv.pl`, etc.).
- This repository does **not** ship the SIRAH force field. You are
  responsible for obtaining it from the SIRAH developers (http://www.sirahff.com)
  and complying with their license and citation requirements.
- `run_setup.sh` checks that `GMXLIB`, if set, exists and that the SIRAH
  folder, `wt416.gro` and `cgconv.pl` are present, and stops with a clear
  error otherwise. The definitive check is always `pdb2gmx` itself at STEP 0.

---

## 7. Analysis

Once production trajectories exist (`4_prd/repN/traj_all_center_PBC.xtc`),
run:

```bash
./run_setup_analysis.sh
```

NOTE: If you used the slurm method, before running the analysis script, 
issue gunzip 001.tpr.gz inside every replica's folder.

See the comments at the top of that script for the individual analyses it
runs and the scripts it calls.

`run_setup_analysis.sh` has its own `EXECUTION_MODE`, independent of
`run_setup.sh`:

- `EXECUTION_MODE=ssh` (default) — dispatches analyses across a pool of
  machines (`MACHINES=(...)`) over SSH, `MAX_JOBS_PER_MACHINE` at a time per
  machine. Requires passwordless SSH to each machine.
- `EXECUTION_MODE=local` — no SSH, no machine pool. Everything runs on the
  current machine instead, `LOCAL_MAX_JOBS` analyses at a time. Set
  `LOCAL_MAX_JOBS=1` to run them strictly one after another, or raise it to
  run a few at once — useful if you only have your own workstation.

Each analysis produces a .png plot alongside its data — per-replica plots land in 4_prd/repN/analysis/<subfolder>/, and pooled analyses (e.g. HDBSCAN, interactions) land in 4_prd/Analysis/<subfolder>/

> **Note on `MACHINES`:** the default hostnames (`bio212`, `bio213`, ...)
> are specific to this repository's original cluster and won't exist
> elsewhere. Replace them with your own reachable hostnames, or just use
> `EXECUTION_MODE=local`.

> **Note on CG analyses:** GROMACS has no default "Backbone" group for CG
> beads (SIRAH 2 backbone beads are `GN` and `GO`), so `index_maker.sh` builds
> it. DSSP is atomistic-only and is not in `ANALYSES[]`; use `sirah_ss` in VMD
> instead (see [section 8](#8-sirah-tools-in-vmd)).

---

## 8. SIRAH Tools in VMD

SIRAH Tools is a set of utilities shipped with the force field
(`sirah.ff/tools/`). Besides the `cgconv.pl` mapping used by `run_setup.sh`,
it includes a VMD plugin, `sirah_vmdtk.tcl`, to visualise SIRAH systems,
calculate secondary structure (`sirah_ss`) and backmap CG to atomistic
(`sirah_backmap`). Type `sirah_help` in the VMD Tcl/Tk console to list all
commands.

Requirements: VMD; AmberTools only for backmapping. Paths below assume you run
from the system root and `0_prep/sirah.ff` is the link created by
`run_setup.sh`.

### 8.1 Check the mapped PDB (before running the simulation)

Load the atomistic input and the CG output side by side:

```bash
vmd -m ab42_1z0q.pdb 0_prep/ab42_CG.pdb
```

Things to check: every residue was mapped, the number of beads matches
`Protein CG beads` in the `run_setup.sh` output, the chains are the ones you
asked for (`N_CHAINS`), and the CG structure overlays the atomistic one.
`cgconv.pl -h` lists its other options.

### 8.2 Load the CG system in VMD

VMD needs a PSF built from the GROMACS topology, and the SIRAH script to set
the right radii, colours and selection macros:

```bash
./0_prep/sirah.ff/tools/g_top2psf.pl -i 0_prep/topol.top -o ab42_CG.psf
vmd ab42_CG.psf 0_prep/ab42_CG.gro -e 0_prep/sirah.ff/tools/sirah_vmdtk.tcl
```

To load a trajectory, the PSF must contain the same beads as the `.xtc`.
`traj_all_center_PBC.xtc` holds only the group set in
`ANALYSIS_CONCAT_CENTER_OUTPUT_GROUP` (default `Protein`), so use a PSF made
from a protein-only topology for it.

### 8.3 Secondary structure

In the VMD Tcl/Tk console, with the trajectory loaded as the top molecule:

```tcl
sirah_ss sel "sirah_protein"
```

This calculates H / E / C for all frames and writes four files:

| File | Content |
|---|---|
| `ss_by_frame.xvg` | H, E, C percentages per frame |
| `ss_by_res.xvg` | H, E, C percentages per residue |
| `ss_global.xvg` | Global H, E, C percentages |
| `ss.mtx` | Per-residue, per-frame assignment matrix |

Useful options:

```tcl
sirah_ss outname {mtx myss.mtx byframe myss_by_frame.xvg}   ;# custom output names
sirah_ss mol 2 load ss.mtx                                   ;# load previous results into molecule 2
```

Once loaded, the Cartoon representation can be coloured by `Structure`
(Graphics → Representations → Coloring Method).

### 8.4 Backmapping (CG → atomistic)

```tcl
sirah_backmap
```

`sirah_backmap` reconstructs pseudo-atomistic coordinates from the CG
structure or from any frame of the trajectory. It calls AmberTools, so
AmberTools must be installed and in your `PATH`. Run `sirah_backmap` without
arguments (or `sirah_help`) to see the available options. The result is a
pseudo-atomistic model whose accuracy is limited by the CG resolution; minimise
it before using it as a starting structure for an all-atom run.

Headless use (no display), e.g. on a cluster node:

```bash
vmd -dispdev text -e script.tcl
```

---

## Citation

If you use SIRAH 2, please cite the SIRAH force-field publications listed at
http://www.sirahff.com. For SIRAH Tools (mapping, backmapping and VMD plugin):
Machado, M. R.; Pantano, S. *Bioinformatics* 2016, 32(10), 1568–1570.
doi: 10.1093/bioinformatics/btw020
