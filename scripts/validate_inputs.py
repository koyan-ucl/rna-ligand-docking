#!/usr/bin/env python3
"""List locally available ligand IDs and test SDF files without running docking."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from docking_project import discover_ligands, find_sdf, validate_sdf

root = Path(__file__).resolve().parents[1]
ligand_dir = root / 'data' / 'ligands'
csv_file = root / 'data' / 'Integrated.csv'
ids = discover_ligands(ligand_dir, csv_file)
print(f'Found {len(ids)} ligand IDs (from CSV or local SDF filenames)')
errors = 0
for ligand_id in ids:
    try:
        path = find_sdf(ligand_id, ligand_dir)
        validate_sdf(path)
        print(f'OK    {ligand_id}: {path.name}')
    except (ValueError, FileNotFoundError) as exc:
        errors += 1
        print(f'ERROR {ligand_id}: {exc}')
print(f'Valid: {len(ids) - errors}; invalid/missing: {errors}')
raise SystemExit(1 if errors or not ids else 0)
