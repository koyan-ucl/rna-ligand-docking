# RNA–Ligand Virtual Screening with BioBB, fpocket and AutoDock Vina

An annotated Jupyter workflow for **RNA binding-pocket identification and virtual screening of previously downloaded small molecules**. Adapted from an Oxford 7LYJ RNA docking notebook. The project prepares an RNA structure, identifies candidate pockets with fpocket, generates a docking box, converts local ligand SDF files to PDBQT, and runs AutoDock Vina through BioBB.

**No online ligand retrieval is performed.** ZINC20 now challenges automated requests with CAPTCHA; the workflow instead accepts your existing PDB/SDF data. This repository does **not** include the user's original molecular files, which must be supplied separately.

> **Research caution:** Docking scores are approximate ranking indicators, **not** experimental binding free energies. In particular, RNA charges/protonation, metal ions, receptor flexibility and Vina's scoring model require careful validation.

## Workflow overview

```text
data/structures/7LYJ.pdb
       |  PDBFixer (prepare RNA)
       v
 results/processed.pdb -----> fpocket -----> filter/select pocket
       |                                      |
       |                                      v
       |                                 results/box.pdb
       v                                      |
 receptor PDBQT                              |
       |                                      |
       +------------------+-------------------+
                          |
data/ligands/*.sdf -> Open Babel (PDB/PDBQT) -> AutoDock Vina
                                               |
                                               v
                           results/batch/results.csv + docked poses
```

## Files

```text
rna-ligand-docking/
├── README.md
├── environment.yml                # Conda packages and required executables
├── .gitignore                     # Excludes docking outputs, not inputs
├── docking_project.py             # SDF discovery, validation, batch workers and scores
├── notebooks/
│   └── RNA_Docking_Workflow.ipynb # Main annotated, clean-output notebook
├── data/
│   ├── Integrated.example.csv     # Example ID list; rename to Integrated.csv
│   ├── structures/                # Your local RNA PDBs (not included)
│   └── ligands/                   # Your local ligand SDFs (not included)
├── scripts/
│   └── validate_inputs.py         # Input-only preflight check
├── tests/
│   └── test_docking_project.py    # Unit tests without docking binaries
└── .github/workflows/test.yml     # Test-only CI; no large docking jobs
```

## 1. Install on Linux (recommended)

```bash
conda env create -f environment.yml
conda activate molecular-jupyter
python -m jupyterlab
```

Or in your **existing** `molecular-jupyter` environment, check that the following executables are on `PATH`:

```bash
which fpocket obabel vina check_structure
vina --version
obabel -V
```

The environment requests Python 3.11, JupyterLab, nglview 3.1.4, OpenMM/PDBFixer, BioBB modules, and the native executables. Actual package solving may vary with platform and time; the environment has **not been solved or run end-to-end** in this build environment. The required executable is named `vina` (the Conda package is also `vina`, **not** `autodock-vina`).

### Important compatibility notes

- With nglview 3.1.4, `pkg_resources` can issue a harmless deprecation warning. Pinning `setuptools<81` retains it.
- BioBB's `str_check_add_hydrogens` `ph` mode was rejected by the underlying `check_structure` installation in the earlier run. The notebook uses **`mode='auto'`** and validates the output PDBQT format. This does **not** prove that RNA partial charges are appropriate. If the output does not contain correctly formatted charges/types, the notebook retains the BioBB file and tries Open Babel as an explicitly labeled fallback; this can change the charge assignment.
- The workflow defaults to serial batch docking. Switching to process workers can help throughput, but verify a small test and per-worker CPU usage first.

## 2. Add your actual molecular files

**Required:** a PDB structure in `data/structures/` and at least one **real** SDF in `data/ligands/`:

```text
data/structures/7LYJ.pdb
data/ligands/ZINC000003830390.sdf
```

Legacy names are also supported:

```text
data/ligands/tempZINC000003830390.sdf
```

Copy the files you already downloaded; the files are **not bundled in this archive**. An SDF containing a ZINC CAPTCHA HTML page is rejected by the preflight check. The project works with local files and never tries to defeat a CAPTCHA.

Optional CSV for a specific order/selection:

```text
# data/Integrated.csv (file content, no comment line)
Column1
ZINC000003830390
ZINC000003831551
```

If `data/Integrated.csv` is absent, the notebook screens the available SDF filenames, excluding the ambiguous generic name `ideal.sdf`. Rename any generic files to their matching ligand IDs. Existing files under an old `sdf/` directory can be copied directly into `data/ligands/`; filenames such as `tempZINC....sdf` are accepted.

Quick validation (no fpocket or Vina required):

```bash
python scripts/validate_inputs.py
```

**Before pushing molecular data to a public GitHub repository:** verify redistributability, collaborator consent, potentially unpublished structures, and GitHub size limits. If the dataset is large or private, keep `data/` on your machine and upload only code/notebook. For reproducibility, document original structure IDs, downloading dates, SDF source and preparation choices.

## 3. Run the notebook

Open `notebooks/RNA_Docking_Workflow.ipynb` in JupyterLab and select the `molecular-jupyter` kernel. Run cells in order. In the configuration cell, change:

- `PDB_CODE` and `PRIMARY_LIGAND_ID` to match your files.
- `SELECTED_POCKET` after inspecting the fpocket JSON/visualization.
- `POCKET_VOLUME_RANGE` and `BOX_OFFSET` if you intend to change the original pocket-selection protocol.
- `MAX_LIGANDS = 5` for a quick smoke test, then `None` for all IDs; start with `MAX_WORKERS = 1`, then try `3`.

Pocket visualization now uses the original **red cavity surfaces and sky-blue contour/mesh surfaces**, not default atom dots. The selected pocket and the docking box are shown as a mesh and connected box edges.

**Receptor preparation:** section 6 first checks BioBB output, prints example invalid PDBQT records, and (only if necessary) retries with Open Babel in rigid-receptor mode (`-xr -xn`). It saves the original BioBB output as `results/prep_receptor_biobb.pdbqt` for comparison. The alternative method is **not chemically equivalent**; check RNA phosphate partial charges and docking sensitivity before publishing scores.

**Outputs:** generated RNA, pocket artifacts and box in `results/`; individual docking in `results/single/`; batch runs and `results/batch/results.csv` in `results/batch/`.

Example `results.csv` columns: `ligand_id`, `status`, `affinity_kcal_mol`, `sdf_path`, `pose_path`, `error`. `ok` means successfully docked, `skipped` means a previously generated result was reused, and `error` gives a human-readable per-ligand failure.

## 4. GitHub and Google Colab

### Publish a new GitHub repository

1. Create an **empty** GitHub repository named `rna-ligand-docking` (do not add a GitHub-generated README, since this project already has one).
2. Unzip this project on your computer, optionally add your downloaded structures and SDFs under `data/`, then run:

```bash
cd rna-ligand-docking
git init
git add .
git commit -m "Initial RNA docking workflow"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/rna-ligand-docking.git
git push -u origin main
```

If you cannot push because you have not authenticated, use GitHub Desktop or set up GitHub authentication first. If desired, add your repository's GitHub Actions status badge after publishing.

### Open from GitHub in Colab

GitHub notebooks can be opened by using Colab's **GitHub** tab or the URL pattern:

`https://colab.research.google.com/github/YOUR_USERNAME/rna-ligand-docking/blob/main/notebooks/RNA_Docking_Workflow.ipynb`

Colab's notebook UI does **not** automatically clone the repository data. In the first code cell set:

```python
GITHUB_REPO_URL = "https://github.com/YOUR_USERNAME/rna-ligand-docking.git"
```

This clones a public repository into `/content/`. If your SDFs are not public/on GitHub, place them in the cloned `data/ligands/` folder (e.g. by uploading from Colab or mounting Drive) before running the notebook.

**Colab environment caveat:** Colab's hosted runtime does not arrive with `fpocket`, `obabel`, `vina` and `check_structure` preinstalled. This repository is packaged for **local Conda first**; on Colab you still need a compatible environment with those native executables and BioBB installed. Simply opening the GitHub notebook in Colab is not enough to execute the docking workflow. No unverified `apt` installation shortcuts are claimed here.

## 5. Scientific and reproducibility notes

The original notebook intentionally removed heterogens and waters. For RNA, removing metal ions such as Mg²⁺ can remove structural stabilizers; decide whether the original method suits your target. Inspect hydrogen addition, phosphate charge assignments and atom types in the PDBQT before production docking. Vina's scoring function is not a substitute for physics-based free-energy calculations, and the original workflow does not fix random seeds. Differences in Vina/fpocket/BioBB versions or pocket selection can change poses and scores.

For an exact re-run, retain the **original downloaded structures**, their provenance, parameters, external executable versions and output logs. The PDBQT helper now checks documented fixed-column charge, atom-type and XYZ fields and prints sample malformed records; it does not establish charge-model validity. CI tests here cover helper code only; they do not validate physical correctness or run docking.

## Software and documentation

- [BioBB VS (Vina and fpocket wrappers)](https://biobb-vs.readthedocs.io/en/latest/)
- [BioBB Chemistry (Open Babel wrapper)](https://biobb-chemistry.readthedocs.io/en/latest/)
- [BioBB Structure Utils](https://biobb-structure-utils.readthedocs.io/en/latest/)
- [fpocket](https://github.com/Discngine/fpocket)
- [AutoDock Vina](https://vina.scripps.edu/)
- [OpenMM and PDBFixer](https://github.com/openmm/pdbfixer)

### License and attribution

This project is adapted from your original notebook. **No repository license is assigned automatically.** Add your preferred license after checking ownership of contributions and any included structures; third-party dependencies retain their own licenses.
