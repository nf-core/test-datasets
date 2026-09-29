# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

# test-datasets: `moleculardynamics`

Test data for the [nf-core/moleculardynamics](https://github.com/nf-core/moleculardynamics) pipeline, which runs molecular dynamics (MD) simulations with GROMACS: structure pre-processing, topology generation, solvation and ionisation, energy minimisation, NVT and NPT equilibration, production MD, PBC removal and RMSD analysis.

> ⚠️ **This is the `moleculardynamics` branch. Do not merge test data into `master`.**

## Content of this repository

```
testdata/
├── experimental_structures/
│   └── 1AKI.pdb                       # Hen egg-white lysozyme crystal structure
├── mdps/
│   ├── em.mdp       em_test.mdp       # Energy minimisation
│   ├── nvt.mdp      nvt_test.mdp      # NVT equilibration
│   ├── npt.mdp      npt_test.mdp      # NPT equilibration
│   └── md.mdp       md_test.mdp       # Production MD
└── samplesheet/
    └── v1.0/
        ├── samplesheet_test.csv       # used by -profile test
        └── samplesheet_test_full.csv  # used by -profile test_full
```

### Structure

`testdata/experimental_structures/1AKI.pdb` is the orthorhombic form of hen egg-white lysozyme, solved by X-ray diffraction at 1.5 Å resolution ([PDB 1AKI](https://www.rcsb.org/structure/1AKI), deposited by Carter et al., 1997).

- It is a single chain (A) of 129 residues with no missing residues or atoms.
- The file also contains 78 crystallographic waters (`HETATM`), which the pipeline removes during pre-processing.
- Lysozyme is the standard system in introductory GROMACS tutorials, which makes it a small, well-understood test case.
- The file is unmodified from the RCSB PDB, which releases its data under CC0 1.0.

### MDP parameter files

Both MDP sets use the same physics:

- non-bonded settings for the CHARMM force field: 1.2 nm cut-offs, force-switched van der Waals, PME electrostatics
- bonds to hydrogen constrained with LINCS
- V-rescale thermostat at 298 K and C-rescale barostat at 1 bar
- position restraints on the protein (`-DPOSRES`) during equilibration

The two sets differ only in run length and output frequency:

| Step | Full (`*.mdp`) | Minimal (`*_test.mdp`) |
|---|---|---|
| Energy minimisation (steepest descent, Fmax < 1000 kJ/mol/nm) | max. 5000 steps | max. 1000 steps |
| NVT equilibration (dt = 2 fs) | 20 ps (10 000 steps) | 1 ps (500 steps) |
| NPT equilibration (dt = 2 fs) | 20 ps (10 000 steps) | 1 ps (500 steps) |
| Production MD (dt = 2 fs) | 50 ps (25 000 steps), frames every 2 ps | 2 ps (1000 steps), frames every 0.2 ps |

The minimal NVT file uses a fixed random seed (`gen_seed = 12345`), so CI runs start from the same velocities every time.

### Samplesheets

Both samplesheets contain one sample, `LYSOZYME_1AKI`, set up with the CHARMM27 force field and a cubic box with a 1.0 nm solute–box distance. They point to the files above using raw GitHub URLs on this branch.

| Samplesheet | Profile | MDP set | Purpose |
|---|---|---|---|
| `samplesheet_test.csv` | `test` | `*_test.mdp` | Fast CI test, run on every pull request (about 1–2 min on 4 CPUs) |
| `samplesheet_test_full.csv` | `test_full` | `*.mdp` | Longer run for the full-size AWS test on releases |

The columns follow the pipeline's input schema:

```
sample,structure,em_mdp,nvt_mdp,npt_mdp,md_mdp,forcefield,box_type,distance_to_box
```

## Usage

```bash
# Minimal test
nextflow run nf-core/moleculardynamics -profile test,docker --outdir <OUTDIR>

# Full-size test
nextflow run nf-core/moleculardynamics -profile test_full,docker --outdir <OUTDIR>
```

To download only this branch:

```bash
git clone https://github.com/nf-core/test-datasets.git --single-branch --branch moleculardynamics
```

## Support

For help, join the `#moleculardynamics` channel on the [nf-core Slack](https://nf-co.re/join/slack). General guidance on test data is in the [nf-core test data guidelines](https://nf-co.re/docs/contributing/test_data_guidelines).
