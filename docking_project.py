"""Local-file utilities for the RNA / ligand BioBB screening notebook.

No network access or ZINC download is performed here. Each worker returns
plain Python values, including a printable error, for safe multiprocessing.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import csv
import re
import shutil
from typing import Any


VINA_DEFAULTS = {"cpu": 1, "exhaustiveness": 8, "num_modes": 9,
                 "min_rmsd": 1, "energy_range": 3}


def find_sdf(ligand_id: str, ligand_dir: Path) -> Path:
    """Locate locally downloaded SDFs with common ZINC naming conventions."""
    ligand_dir = Path(ligand_dir)
    ligand_id = str(ligand_id).strip()
    if not ligand_id or Path(ligand_id).name != ligand_id:
        raise ValueError(f"Invalid ligand ID: {ligand_id!r}")
    for name in (f"{ligand_id}.sdf", f"temp{ligand_id}.sdf",
                 f"{ligand_id.lower()}.sdf", f"temp{ligand_id.lower()}.sdf"):
        path = ligand_dir / name
        if path.is_file():
            return path
    raise FileNotFoundError(
        f"SDF for {ligand_id} not found in {ligand_dir}. "
        f"Expected {ligand_id}.sdf or temp{ligand_id}.sdf"
    )


def discover_ligands(ligand_dir: Path, csv_file: Path | None = None,
                     csv_column: str = "Column1") -> list[str]:
    """Use IDs from Integrated.csv if provided, otherwise discover SDFs."""
    ligand_dir = Path(ligand_dir)
    if csv_file is not None and Path(csv_file).is_file():
        with Path(csv_file).open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            if csv_column not in (reader.fieldnames or []):
                raise ValueError(f"CSV must contain column {csv_column!r}; "
                                 f"found {reader.fieldnames}")
            ids = [str(row[csv_column]).strip() for row in reader
                   if row.get(csv_column) and str(row[csv_column]).strip()]
    else:
        ids = []
        for path in sorted(ligand_dir.glob("*.sdf")):
            name = path.stem
            if name.lower().startswith("temp"):
                name = name[4:]
            if name.lower() == "ideal":
                continue  # generic filename does not encode a ligand identity
            ids.append(name)
        ids.sort()  # consistent ID ordering regardless of "temp" prefix
    return list(dict.fromkeys(ids))


def validate_sdf(path: Path) -> None:
    """Catch HTML CAPTCHA pages and obvious invalid SDFs early."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty SDF: {path}")
    raw = path.read_bytes()
    if b"<html" in raw[:4096].lower() or b"verify you're human" in raw.lower():
        raise ValueError(f"{path} is HTML (probably a CAPTCHA page), not SDF")
    if b"M  END" not in raw:
        raise ValueError(f"{path} does not contain an SDF/MOL 'M  END' record")
    if raw.count(b"M  END") > 1:
        raise ValueError(f"{path} contains multiple molecules; split them first")


def validate_pdbqt(path: Path) -> dict[str, int]:
    """Check that the receptor file has actual PDBQT atom records and charges.

    This is a format sanity check, NOT a chemical validation of RNA charges.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    lines = [line for line in path.read_text(errors="replace").splitlines()
             if line.startswith(("ATOM  ", "HETATM"))]
    if not lines:
        raise ValueError(f"{path}: no ATOM/HETATM records")
    charges = 0
    for line in lines:
        if len(line) >= 79:
            try:
                float(line[70:76]); atom_type = line[77:79].strip()
                if atom_type:
                    charges += 1
            except ValueError:
                pass
    if charges != len(lines):
        raise ValueError(
            f"{path}: {charges}/{len(lines)} atoms have readable PDBQT "
            "charges and atom types. Recheck RNA receptor preparation."
        )
    return {"atom_count": len(lines), "charged_atom_count": charges}


def best_affinity(pdbqt: Path) -> float | None:
    """Return best Vina score (kcal/mol) from PDBQT REMARK VINA RESULT."""
    path = Path(pdbqt)
    if not path.exists():
        return None
    match = re.search(
        r"^REMARK\s+VINA\s+RESULT:\s*([-+]?\d+(?:\.\d+)?)",
        path.read_text(errors="replace"), flags=re.MULTILINE,
    )
    return float(match.group(1)) if match else None


def _dock_worker(job: tuple[str, str, str, str, str, dict[str, Any], bool]) -> dict:
    """Module-level subprocess-friendly worker. Never returns raw exceptions."""
    ligand_id, sdf, receptor, box, out_dir, vina_props, skip_existing = job
    item_dir = Path(out_dir) / ligand_id
    item_dir.mkdir(parents=True, exist_ok=True)
    pdb = item_dir / "ligand.pdb"
    pdbqt = item_dir / "ligand.pdbqt"
    out_pdbqt = item_dir / "docked.pdbqt"
    log = item_dir / "vina.log"
    result = {"ligand_id": ligand_id, "status": "error", "affinity_kcal_mol": None,
              "sdf_path": sdf, "pose_path": str(out_pdbqt), "error": ""}
    try:
        if skip_existing:
            affinity = best_affinity(out_pdbqt)
            if affinity is not None:
                return {**result, "status": "skipped", "affinity_kcal_mol": affinity}
        validate_sdf(Path(sdf))
        from biobb_chemistry.babelm.babel_convert import babel_convert
        from biobb_vs.vina.autodock_vina_run import autodock_vina_run
        if babel_convert(input_path=sdf, output_path=str(pdb), properties={
            "input_format": "sdf", "output_format": "pdb", "binary_path": "obabel"
        }) != 0:
            raise RuntimeError("Open Babel SDF→PDB conversion failed")
        if babel_convert(input_path=str(pdb), output_path=str(pdbqt), properties={
            "input_format": "pdb", "output_format": "pdbqt", "binary_path": "obabel"
        }) != 0:
            raise RuntimeError("Open Babel PDB→PDBQT conversion failed")
        rc = autodock_vina_run(
            input_ligand_pdbqt_path=str(pdbqt),
            input_receptor_pdbqt_path=receptor,
            input_box_path=box,
            output_pdbqt_path=str(out_pdbqt),
            output_log_path=str(log),
            properties=dict(vina_props),
        )
        if rc != 0:
            raise RuntimeError(f"Vina exited with status {rc}")
        affinity = best_affinity(out_pdbqt)
        if affinity is None:
            raise RuntimeError("No Vina score found in output pose; check vina.log")
        return {**result, "status": "ok", "affinity_kcal_mol": affinity}
    except Exception as exc:
        # Return printable values rather than potentially unpicklable BioBB errors.
        return {**result, "error": f"{type(exc).__name__}: {exc}"}


def run_batch(ligand_ids: list[str], ligand_dir: Path, receptor: Path,
              box: Path, out_dir: Path, max_workers: int = 1,
              vina_props: dict | None = None, skip_existing: bool = True) -> list[dict]:
    """Dock local SDF files; errors are recorded per ligand, not lost in Pool.map."""
    if max_workers < 1:
        raise ValueError("max_workers must be >= 1")
    validate_pdbqt(Path(receptor))
    if not Path(box).is_file():
        raise FileNotFoundError(box)
    jobs = []
    records = []
    for ligand_id in ligand_ids:
        try:
            sdf = find_sdf(ligand_id, ligand_dir)
            jobs.append((ligand_id, str(sdf.resolve()), str(Path(receptor).resolve()),
                         str(Path(box).resolve()), str(Path(out_dir).resolve()),
                         vina_props or VINA_DEFAULTS, skip_existing))
        except Exception as exc:
            records.append({"ligand_id": ligand_id, "status": "error",
                            "affinity_kcal_mol": None, "sdf_path": "", "pose_path": "",
                            "error": f"{type(exc).__name__}: {exc}"})
    if max_workers == 1:
        records.extend(_dock_worker(job) for job in jobs)
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as pool:
            futures = [pool.submit(_dock_worker, job) for job in jobs]
            for future in as_completed(futures):
                records.append(future.result())
    records.sort(key=lambda row: (row["affinity_kcal_mol"] is None,
                                   row["affinity_kcal_mol"] if row["affinity_kcal_mol"] is not None else 0,
                                   row["ligand_id"]))
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    with (Path(out_dir) / "results.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=(
            "ligand_id", "status", "affinity_kcal_mol", "sdf_path", "pose_path", "error"))
        writer.writeheader()
        writer.writerows(records)
    return records
