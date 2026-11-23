# SIRAH 2 — Technical and Theoretical Notes

Companion notes to the README. The README says *how to run* the workflow; these
notes explain *what each step does and why*, so you can check results, debug
problems and adapt the workflow to another system.

> Items marked **(verify)** come from general knowledge of SIRAH and should be
> confirmed against the SIRAH documentation (http://www.sirahff.com) or your
> own tests before being relied on.

---

## 1. What coarse-graining is

An all-atom (AA) simulation of Aβ42 in water is dominated by solvent: a
42-residue peptide has ~600 atoms, while the box holds tens of thousands of
water atoms. Coarse-graining (CG) groups several atoms into one **bead**, which
reduces the number of particles and removes the fastest motions. That allows:

- a larger time step (0.02 ps here vs 0.002 ps in the AA workflow);
- fewer interactions per step;
- longer sampled times for the same compute.

The price is lost detail: bonded and non-bonded terms are *effective*
potentials tuned to reproduce selected properties, so a CG model is only as good
as its parametrisation and the questions you ask of it.

**Rule of thumb:** CG is good for conformational ensembles, collective
behaviour and association; it is not for questions that depend on precise atomic
geometry (exact hydrogen-bond patterns, side-chain rotamers, ligand contacts).

---

## 2. The SIRAH 2 force field

SIRAH (Southamerican Initiative for a Rapid and Accurate Hamiltonian) is a
residue-based CG force field for proteins, nucleic acids, lipids, solvent and
ions. Points relevant to this workflow:

### 2.1 Protein representation

- Beads are **mapped from atom positions** (not from centres of mass of
  arbitrary groups): each bead sits on, or is built from, a defined set of heavy
  atoms. This is what makes mapping *and* backmapping possible with simple
  geometry, and is why RMSD, Rg, etc. can be measured directly on CG
  trajectories.
- The backbone keeps explicit beads for the amide nitrogen and carbonyl oxygen
  (`GN`, `GO`) together with a Cα-type bead. The `Backbone` index group used in
  the analysis is built from `GN` and `GO` for this reason.
- Side chains use a residue-dependent number of beads (small residues few,
  aromatic/charged residues more).
- Partial charges are explicit on the backbone and side-chain beads. Secondary
  structure arises through these electrostatic interactions together with the
  bonded terms, rather than through a dedicated hydrogen-bond potential
  **(verify)**.

### 2.2 Solvent: WT4

SIRAH uses a four-bead water model, **WT4**: one WT4 "molecule" is a small
tetrahedral cluster that represents roughly **11 real water molecules**
**(verify)**, with two positive and two negative beads (`WN1`, `WN2`, `WP1`,
`WP2`), so it carries a dipole-like charge distribution and screens
electrostatics explicitly.

Consequences you see in the scripts:

- `solvate` is fed with `sirah.ff/wt416.gro` (a pre-equilibrated WT4 box).
- The number of WT4 molecules is counted from one atom type (`WP1`) in the
  `.gro` file. `set_molecule_count` rewrites `[ molecules ]` in `topol.top` so it
  matches the `.gro`.
- Water closer than `CLOSE_WATER_CUTOFF` (0.3 nm) to the protein is removed:
  WT4 is bigger than a water molecule, so insertion with the standard solvate
  routine leaves steric clashes next to the solute.

### 2.3 Ions

Ions are **hydrated**: `NaW` and `ClW` represent Na⁺ and Cl⁻ together with their
first hydration shell. They replace WT4 molecules in `genion`
(`ION_SOLVENT_GROUP=WT4`). Charge neutralisation and the salt concentration
(`ION_CONC`) work as in AA, but the counts are of hydrated ions.

### 2.4 Time step and time scale

The standard time step is 20 fs (`NPT_DT` / `PROD_DT` = 0.02 ps). CG
dynamics are smoother than atomistic ones, so simulated time is **not**
equivalent to real time; do not directly compare CG and AA times without
checking the literature **(verify)**.

---

## 3. What each stage of `run_setup.sh` does

| Step | Tool | Purpose |
|---|---|---|
| 0 | `cgconv.pl`, `pdb2gmx` | Map atoms → beads; build the topology. `-ter` answers (`TERMINAL_N/C`) set the termini per chain |
| 1 | `editconf` | Place the solute in a box; `BOX_PADDING` is the minimum solute–edge distance. The box must be larger than twice the interaction cut-off plus any extended conformation, or the peptide interacts with its own image |
| 2 | `solvate` + `select` | Fill with WT4, then drop WT4 within 0.3 nm of the protein |
| 3 | `genion` | Neutralise and add salt (hydrated ions) |
| 3b | `make_ndx`, `genrestr` | Build the index and the backbone position-restraint files (`bkbres*.itp`, `bkbres_soft*.itp`) |
| 4 | `min_sirah2.sh` | Two-stage energy minimisation (`em1` then `em2`) to remove clashes |
| 5 | `eq_NVT_sirah2.sh` | Equilibrate temperature at fixed volume, with restraints |
| 6 | NPT (`eq2_CGPROT.mdp`) | Equilibrate pressure/density; restraints optional (`POSRES_NPT`) |
| 7 | Production (`md_CGPROT.mdp`) | Unrestrained simulation, in `PROD_SEGMENTS` blocks, `REPLICATES` times |

### Why restrain the backbone during equilibration?

After minimisation, solvent is not yet relaxed around the solute. Restraining
the backbone beads lets water and ions equilibrate without the peptide
distorting. `bkbres` is the stronger restraint and `bkbres_soft` the weaker one
(used in the later stage). For multi-chain systems one pair is generated per
chain, because `topol.top` includes one `.itp` per chain.

### Why segments and replicates?

- **Segments:** each block runs `nsteps × dt` and continues from the previous
  one. Segment *n* starts at `tinit = nsteps × dt × (n − 1)`. This fits cluster
  wall-time limits and makes runs restartable.
- **Replicates:** IDPs such as Aβ42 do not have one native state. A single
  trajectory can get trapped, so several independent runs from the same NPT
  endpoint (different velocities) are needed to estimate uncertainty.

---

## 4. Mapping: from atoms to beads (`cgconv.pl`)

`cgconv.pl -i atomistic.pdb -o cg.pdb` reads the heavy atoms and writes the
bead positions. Since beads are defined from atom positions, the CG PDB can be
overlaid on the atomistic structure for a sanity check (see section 5.2).

What can go wrong:

- **Missing atoms** in the input PDB → missing or misplaced beads.
- **Non-standard residues or altered naming** → residue not recognised.
- **Several chains / models (NMR files):** NMR PDBs contain many models. Use one
  model, or the CG PDB will contain all of them. `N_CHAINS` keeps only the first
  N chains from the output.
- **Protonation:** mapping uses the protonation state encoded in the input
  (e.g. residue names such as HIS variants). Check histidines.

Always confirm that the number of beads printed by `run_setup.sh`
(`Protein CG beads`) is what you expect for your sequence.

---

## 5. Working with SIRAH Tools in VMD

SIRAH Tools = `cgconv.pl` (mapping), `g_top2psf.pl` (topology → PSF) and the VMD
plugin `sirah_vmdtk.tcl` (visualisation, analysis, backmapping). Reference:
Machado & Pantano, *Bioinformatics* 2016, 32, 1568–1570.

### 5.1 Why a plugin is needed

VMD does not know SIRAH bead types: it would assign default radii and colours,
and it would not recognise SIRAH residues as protein. The plugin sets the van
der Waals radii, residue types and selection macros (e.g. `sirah_protein`), so
that cartoon-style representations, selections and colouring work. In GROMACS,
`residuetypes.dat` and `specbond.dat` play the same role, which is why
`run_setup.sh` links them into every stage folder.

### 5.2 Checking the mapping

```bash
vmd -m ab42_1z0q.pdb 0_prep/ab42_CG.pdb
```

Checklist:

1. Both molecules overlay (beads sit on the atoms they were built from).
2. No residue is missing or doubled.
3. Chain count = `N_CHAINS`.
4. Bead count matches the log.
5. Termini look right (charged vs neutral choice in `TERMINAL_N/C`).

### 5.3 Loading a system or trajectory

```bash
./0_prep/sirah.ff/tools/g_top2psf.pl -i 0_prep/topol.top -o ab42_CG.psf
vmd ab42_CG.psf 0_prep/ab42_CG.gro -e 0_prep/sirah.ff/tools/sirah_vmdtk.tcl
```

The PSF holds the **topology** (atoms, bonds); the `.gro`/`.xtc` hold the
**coordinates**. They must have the same number of particles in the same order.
`traj_all_center_PBC.xtc` contains only the group chosen in
`ANALYSIS_CONCAT_CENTER_OUTPUT_GROUP` (default `Protein`), so pair it with a
protein-only PSF **(verify)**. A mismatch usually appears as VMD refusing the
frames or showing scrambled bonds.

### 5.4 Secondary structure with `sirah_ss`

Atomistic DSSP needs backbone N, Cα, C, O and hydrogens; SIRAH lacks them, so
DSSP cannot be used. `sirah_ss` assigns **H** (helix), **E** (strand) and **C**
(coil) from the geometry of the CG backbone. The SIRAH developers report high
agreement with STRIDE on atomistic references.

```tcl
sirah_ss sel "sirah_protein"
```

Outputs:

| File | Use |
|---|---|
| `ss_by_frame.xvg` | % H/E/C versus time — shows folding/unfolding events |
| `ss_by_res.xvg` | % H/E/C per residue — shows which residues form β-sheet |
| `ss_global.xvg` | Overall H/E/C content |
| `ss.mtx` | Full matrix (residue × frame); reload with `sirah_ss mol N load ss.mtx` |

How to read them for an IDP like Aβ42:

- Expect mostly C with transient E; a persistent E in particular segments is
  what you look for when studying β-sheet formation.
- Plot `ss_by_res` per replicate and compare: agreement between replicates is
  a sign of convergence.
- A three-state assignment (H/E/C) is coarser than DSSP's eight states; do not
  compare percentages with DSSP numbers without caveats.

### 5.5 Backmapping: from CG to atomistic

**Idea.** Because every bead is defined from specific atoms, the reverse
operation is a geometric reconstruction: the missing atoms are placed using
the bead positions and ideal internal geometry, giving a **pseudo-atomistic**
model. In SIRAH Tools this is `sirah_backmap`, which relies on AmberTools to
complete/optimise the structure.

```tcl
sirah_backmap
```

Points to remember:

- The model is **only as good as the CG information**: side-chain rotamers and
  hydrogen-bond details are not recovered from the CG trajectory; they are
  rebuilt.
- Always **energy-minimise** the result (and ideally equilibrate it) before using
  it as an all-atom starting structure; clashes and strained angles are expected.
- You can backmap any frame, so a common use is to pick representative
  structures (e.g. cluster centres from the HDBSCAN analysis) and backmap them
  to run all-atom simulations or atomistic analyses on them (multiscale
  workflow).
- For the full list of options, use `sirah_help` or run the command without
  arguments in the VMD Tcl/Tk console.

---

## 6. Analysis in CG: what changes

| Analysis | In CG |
|---|---|
| Rg | Works directly; beads are mapped from atoms, but the CG Rg is not identical to the AA one — compare CG vs CG |
| SASA | Depends on bead radii; absolute values differ from AA. Use for relative comparisons |
| Secondary structure | `sirah_ss` in VMD instead of DSSP |
| Contacts / interactions | Defined between beads; cut-offs need to be adapted to bead size |
| Clustering (HDBSCAN) | Works on bead coordinates; remember `nowater` files are needed |
| Index groups | No default `Backbone` group; `GN`+`GO` are built in `index_maker.sh` from `system_groups.cfg` |

**Do not mix CG and AA numbers in the same plot without stating it.** Compare
trends, not absolute values.

---

## 7. Practical checks and troubleshooting

| Symptom | Likely cause |
|---|---|
| `pdb2gmx` fails at STEP 0 | Wrong `FF_NAME`/`WM`, broken `sirah.ff` link, or `TERMINAL_N/C` menu numbers invalid for your system |
| Requested N chains but fewer found | Chain IDs missing/repeated in the PDB; check `TER` records |
| `grompp` warns about charge/molecule counts | `[ molecules ]` does not match the `.gro` (check WT4, `NaW`, `ClW` counts) |
| Minimisation does not converge | Clashes between protein and WT4; increase `CLOSE_WATER_CUTOFF` slightly or check the mapping |
| No `bkbres*.itp` | Index/restraint step failed — check `make_index_restraints_sirah2.sh` output and chain naming |
| VMD shows wrong bonds/colours | PSF does not match the coordinates, or `sirah_vmdtk.tcl` not loaded |
| Trajectory looks fragmented | Periodic boundaries; use the concatenated, centred `traj_all_center_PBC.xtc` |

---

## 8. Glossary

- **Bead** — CG particle representing a group of atoms.
- **WT4** — SIRAH four-bead water model (one molecule = a cluster of waters).
- **NaW / ClW** — hydrated sodium / chloride ions.
- **GN / GO** — backbone amide N and carbonyl O beads.
- **Mapping** — atoms → beads. **Backmapping** — beads → atoms.
- **Pseudo-atomistic** — atomistic model rebuilt from CG coordinates, not
  from an atomistic simulation.
- **PSF** — VMD topology file (atoms, bonds, types).
- **Segment** — one continuation block of `nsteps × dt`.
- **Replicate** — independent production run from the same equilibrated state.

---

## References

- Machado, M. R.; Barrera, E. E.; Klein, F.; Sóñora, M.; Silva, S.; Pantano, S.
  The SIRAH 2.0 Force Field: Altius, Fortius, Citius. *J. Chem. Theory Comput.*
  2019, 15, 2719–2733.
- Machado, M. R.; Pantano, S. SIRAH Tools: mapping, backmapping and
  visualization of coarse-grained models. *Bioinformatics* 2016, 32, 1568–1570.
- SIRAH documentation and tutorials: http://www.sirahff.com
