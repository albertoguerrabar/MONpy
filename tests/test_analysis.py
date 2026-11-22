"""Regression checks for the published acrolein example."""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "acrolein"


class TestAcrolein(unittest.TestCase):
    def run_case(self, config_text):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            config_dir = tmp / "calculation"
            config_dir.mkdir()
            shutil.copyfile(EXAMPLE / "acro.ndo", config_dir / "acro.ndo")
            config = config_dir / "monpy.toml"
            config.write_text(config_text, encoding="utf-8")
            proc = subprocess.run(
                [sys.executable, "-B", str(ROOT / "main.py"), "--config", str(config)],
                cwd=tmp, stdin=subprocess.DEVNULL, capture_output=True, text=True,
                env={**os.environ, "MPLCONFIGDIR": str(tmp / "mpl"),
                     "XDG_CACHE_HOME": str(tmp / "cache")},
                check=False,
            )
            results = config_dir / "results"
            if not results.is_dir():
                return proc, [], set()
            with (results / "estados.csv").open(encoding="utf-8", newline="") as stream:
                rows = list(csv.DictReader(stream))
            return proc, rows, {path.name for path in results.iterdir()}

    def test_reference_results(self):
        proc, rows, files = self.run_case((EXAMPLE / "monpy.toml").read_text())
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual([(r["spin"], r["state"]) for r in rows], [("singlet", "1"), ("singlet", "2")])
        for row, expected in zip(rows, [0.081304, 0.060786]):
            self.assertEqual(row["status"], "ok")
            self.assertAlmostEqual(float(row["lambda_CT"]), expected, places=6)
            self.assertAlmostEqual(
                float(row["P[1,2]"]) + float(row["P[2,1]"]),
                float(row["lambda_CT"]),
            )
        self.assertEqual(files, {
            "monpy_results.txt", "estados.csv", "CT_singlet_1.png", "CT_singlet_2.png",
        })

    def test_rejects_overlapping_fragments(self):
        config = (EXAMPLE / "monpy.toml").read_text().replace('B = "5-8"', 'B = "4-8"')
        proc, _, files = self.run_case(config)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Átomo repetido en la partición: 4", proc.stderr)
        self.assertFalse(files)


if __name__ == "__main__":
    unittest.main()
