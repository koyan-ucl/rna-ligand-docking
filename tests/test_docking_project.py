from pathlib import Path
import pytest
from docking_project import find_sdf, discover_ligands, validate_sdf, best_affinity


def test_find_sdf_supports_legacy_temp_name(tmp_path):
    (tmp_path / "tempZINC123.sdf").write_text("test\nM  END\n$$$$\n")
    assert find_sdf("ZINC123", tmp_path).name == "tempZINC123.sdf"


def test_discover_ligands_finds_both_names(tmp_path):
    for filename in ("tempZINC001.sdf", "ZINC002.sdf", "ideal.sdf"):
        (tmp_path / filename).write_text("M  END\n")
    assert discover_ligands(tmp_path) == ["ZINC001", "ZINC002"]


def test_validate_sdf_rejects_captcha_html(tmp_path):
    path = tmp_path / "ZINC001.sdf"
    path.write_text("<html><title>ZINC - Verify you're human</title></html>")
    with pytest.raises(ValueError, match="HTML"):
        validate_sdf(path)


def test_best_affinity(tmp_path):
    path = tmp_path / "docked.pdbqt"
    path.write_text("MODEL 1\nREMARK VINA RESULT: -7.600 0.000 0.000\nENDMDL\n")
    assert best_affinity(path) == -7.6
